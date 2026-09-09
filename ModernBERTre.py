import os
import random, numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from transformers import AutoModel, AutoConfig, Trainer, TrainingArguments,AutoModelForSequenceClassification
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, average_precision_score, matthews_corrcoef
import torch.nn as nn
from scipy.special import softmax
import argparse
import torch.nn.functional as F
# 设置环境变量
os.environ['TOKENIZERS_PARALLELISM'] = 'false'  # 防止 tokenizer 多进程警告
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
torch.backends.cudnn.benchmark = True
# output_dir = "./checkpoints_five"
# os.makedirs(output_dir, exist_ok=True)
def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)
    # 禁用TF32确保FP32精确计算
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
set_seed(42)   # 文章中注明此值
# 加载数值型向量数据

# 构造数据集
class VecDataset(Dataset):
    def __init__(self, feats, labels):
        self.feats = torch.tensor(feats, dtype=torch.float32)
        self.labels = torch.tensor(labels, dtype=torch.long)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        # return {"vec": self.feats[idx], "labels": self.labels[idx]}
        return {
            "vec": self.feats[idx],
            "labels": self.labels[idx]
        }

class ImprovedGODown(nn.Module):
    def __init__(self, go_dim, hidden_dims, output_dim):
        super().__init__()
        layers = []
        prev_dim = go_dim

        # 构建渐进降维层
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.GELU(),
                nn.Dropout(0.2)
            ])
            prev_dim = hidden_dim

        # 输出层
        layers.append(nn.Linear(prev_dim, output_dim))

        self.mlp = nn.Sequential(*layers)

    def forward(self, x):
        return self.mlp(x)
# 修改 ModernBERT 模型以接受数值型向量数据
class ModernBERTForVectorClassification(torch.nn.Module):
    def __init__(self, model_path, num_labels):
        super(ModernBERTForVectorClassification, self).__init__()
        # self.bert = AutoModel.from_pretrained(model_path)
        self.bert = AutoModelForSequenceClassification.from_pretrained(model_path, num_labels=num_labels,attn_implementation="eager")
        # self.dropout = torch.nn.Dropout(0.1)
        # self.classifier = torch.nn.Linear(self.bert.config.hidden_size, num_labels)
        self.go_down = ImprovedGODown(go_dim=13805, hidden_dims=[2048, 512], output_dim=64 * 768)
        self.expr_proj = nn.Linear(32, 32 * 768)  # 关键修改
        # self.input_projection = torch.nn.Linear(192, self.bert.config.hidden_size)
        self.pos_embed = nn.Parameter(torch.randn(1, 192, 768) * 0.02)  # 32+64=96×2=192

    def forward(self, vec, labels=None):
        batch_size, device = vec.size(0), vec.device
    # def forward(self, vec, labels=None):
        e_src = vec[:, :32]  # 源基因表达特征
        g_src = vec[:, 32:32 + 13805]  # 源GO特征
        e_tgt = vec[:, 32 + 13805:32 + 13805 + 32]  # 目标基因表达特征
        g_tgt = vec[:, 32 + 13805 + 32:]  # 目标GO特征
        src_g = self.go_down(g_src).view(batch_size, 64, -1)
        src_e = self.expr_proj(e_src).view(batch_size, 32, -1)
        tgt_e = self.expr_proj(e_tgt).view(batch_size, 32, -1)
        tgt_g = self.go_down(g_tgt).view(batch_size, 64, -1)
        src = torch.cat([src_e, src_g], dim=1)  # [B, 96, d_model]
        tgt = torch.cat([tgt_e, tgt_g], dim=1)  # [B, 96, d_model]
        seq = torch.cat([src, tgt], dim=1)  # [B, 192, d_model]
        seq = seq + self.pos_embed.to(device)
        attention_mask = torch.ones(batch_size, 192, device=device)
        # inputs = self.input_projection(vec.unsqueeze(1))
        # 将数值型向量数据通过一个线性层映射到 BERT 的输入维度
        # outputs = self.bert(inputs_embeds=inputs)
        outputs = self.bert(inputs_embeds=seq, labels=labels, attention_mask=attention_mask)

        # last_hidden_state = outputs.last_hidden_state  # 获取最后一层隐藏状态
        # cls_token = last_hidden_state[:, 0, :]  # 提取 [CLS] 标记（索引0）
        #
        # pooled_output = self.dropout(cls_token)  # 应用Dropout
        # logits = self.classifier(pooled_output)
        #
        # loss = None
        # if labels is not None:
        #     loss = torch.nn.CrossEntropyLoss()(logits, labels)
        # return {"loss": loss, "logits": logits}
        return outputs





# 定义评估指标
def compute_metrics(eval_pred):
    logits, labels = eval_pred
    preds = np.argmax(logits, axis=-1)
    probs = softmax(logits, axis=-1)[:, 1]
    acc = accuracy_score(labels, preds)
    precision_macro = precision_score(labels, preds, average='macro')
    neg_recall = recall_score(labels, preds, pos_label=0)
    pos_recall = recall_score(labels, preds, pos_label=1)
    g_mean = (neg_recall * pos_recall) ** (1 / 2)
    f1 = f1_score(labels, preds, average='macro')
    auroc = roc_auc_score(labels, probs)
    aupr = average_precision_score(labels, probs)
    mcc = matthews_corrcoef(labels, preds)
    return {
        "A_accuracy": acc,
        "B_precision": precision_macro,
        "C-G-mean": g_mean,
        "D_f1": f1,
        "E_mcc": mcc,
        "F_auroc": auroc,
        "G_aupr": aupr
    }
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--txt_input", required=True, help="nature_amps.csv")
    parser.add_argument("--array_input", required=True, help="Artifical_amps.csv")
    parser.add_argument("--output_dir", required=True, help="Artifical_amps.csv")
    parser.add_argument("--count_fold", required=True, help="Artifical_amps.csv")
    # parser.add_argument("--output_dir", required=True, help="Artifical_amps.csv")
    args = parser.parse_args()
    expr_go = np.load(args.array_input).astype(np.float32)
    # expr_matrix = expr_go[:, :32]
    # go_matrix = expr_go[:, 32:]
    # expr_go = np.concatenate((expr_matrix, go_matrix), axis=1)
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    X, y = [], []
    val_idx = []  # 记录验证集样本下标
    with open(args.txt_input) as f:
        for idx, line in enumerate(f):
            parts = line.strip().split('\t')
            src_idx = int(parts[2])
            tgt_idx = int(parts[3])
            label = 1 if parts[4] == '"regulation"' else 0
            sample_tp = parts[5]

            vec = np.concatenate([expr_go[src_idx], expr_go[tgt_idx]])
            X.append(vec.astype(np.float32))
            y.append(label)

            # 划分验证集
            if sample_tp == args.count_fold:
                val_idx.append(idx)

    X = np.stack(X)
    y = np.array(y, dtype=np.int64)

    # 划分训练集和验证集
    val_idx = np.array(val_idx)
    train_mask = np.ones(len(y), dtype=bool)
    train_mask[val_idx] = False
    X_train, y_train = X[train_mask], y[train_mask]
    X_val, y_val = X[val_idx], y[val_idx]
    train_dataset = VecDataset(X_train, y_train)
    val_dataset = VecDataset(X_val, y_val)
    model_path = "/home/lthpc/deepseek/ModernBERT-base"
    model = ModernBERTForVectorClassification(model_path, num_labels=2)
    # 定义训练参数
    training_args = TrainingArguments(
        output_dir=output_dir,
        num_train_epochs=30,
        per_device_train_batch_size=16,
        per_device_eval_batch_size=32,
        learning_rate=1e-6,
        weight_decay=0.01,
        lr_scheduler_type="cosine",
        warmup_ratio=0.1,
        fp16=torch.cuda.is_available(),
        evaluation_strategy="epoch",
        save_strategy="epoch",
        save_total_limit=2,
        load_best_model_at_end=True,
        metric_for_best_model="A_accuracy",
        report_to="none",
        seed=42
    )
    # 初始化 Trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        compute_metrics=compute_metrics
    )

    # 训练模型
    trainer.train()

    # 保存最优模型
    best_dir = os.path.join(output_dir, "best")
    trainer.save_model(best_dir)
    print(f"训练完成，最佳模型已保存到 {best_dir}")