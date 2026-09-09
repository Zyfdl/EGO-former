
import random, numpy as np
import os
import argparse
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
# ★ 必须放在 import transformers 之前
os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"   # 去掉 TF 无关警告
# torch.set_float32_matmul_precision('high')  # 开启 TF32
from transformers import AutoConfig, Trainer, TrainingArguments
from sklearn.metrics import (accuracy_score, recall_score, f1_score,
                             confusion_matrix, roc_auc_score,
                             matthews_corrcoef, average_precision_score,precision_score,)
import torch
import torch.nn as nn
import torch.nn.functional as F

from scipy.special import softmax

torch.backends.cudnn.benchmark = True

from torch.nn import TransformerEncoder, TransformerEncoderLayer
# ---------- 1. 读取原始矩阵 ----------

# expr_matrix = expr_go[:, :32]
# go_matrix = expr_go[:, 32:]
# expr_matrix = np.load('../../imagepy/express_logone.npy').astype(np.float32)
# expr_go = np.concatenate((expr_matrix,go_matrix), axis=1)
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

class VecDataset(torch.utils.data.Dataset):
        def __init__(self, feats, labels):
            self.feats  = torch.tensor(feats, dtype=torch.float32)
            self.labels = torch.tensor(labels, dtype=torch.long)
        def __len__(self): return len(self.labels)
        def __getitem__(self, idx):
            return {"vec": self.feats[idx], "labels": self.labels[idx]}


class GODown(nn.Module):
    def __init__(self, go_dim=13805, hidden_dims=[2048, 512], output_dim=64 * 128):
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

class EGO_former(nn.Module):
    def __init__(self,
                 expr_dim=32,
                 go_dim=13805,
                 d_model=128,  # 特征维度
                 nhead=8,  # 注意力头数
                 num_layers=4,  # Transformer层数
                 max_go_tokens=64,  # GO最大token数
                 dropout=0.15,  # Dropout率
                 num_labels=2,
                 go_hidden_dims=[2048, 512]):  # 分类标签数
        super().__init__()

        # 1. GO降维层（13805维→64个token）
        # self.go_down = nn.Linear(go_dim, max_go_tokens * d_model)
        self.go_down = GODown(go_dim, go_hidden_dims, max_go_tokens * d_model)

        # 2. 基因表达投影层（32维→32个独立token）
        self.expr_proj = nn.Linear(expr_dim, 32 * d_model)  # 关键修改

        # 3. 位置编码（适应新token数量）
        self.pos_embed = nn.Parameter(torch.randn(1, 192, d_model) * 0.02 )# 32+64=96×2=192

        # 4. Transformer编码器
        encoder_layer = TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            activation='gelu',
            batch_first=True
        )
        self.encoder = TransformerEncoder(encoder_layer, num_layers)

        # # 5. 门控特征融合层
        # self.fusion_gate = nn.Sequential(
        #     nn.Linear(d_model * 2, d_model),  # 输入: 双倍维度
        #     nn.ReLU(),
        #     nn.Dropout(dropout),
        #     nn.Linear(d_model, 1),  # 输出单值门控权重
        #     nn.Sigmoid()
        # )
        # self.bilinear_proj_u = nn.Linear(d_model, 96)  # src_embed -> k
        # self.bilinear_proj_v = nn.Linear(d_model, 96)  # tgt_embed -> k
        # 6. 分类头
        self.classifier = nn.Sequential(
            nn.LayerNorm(3*d_model),  # 改成 2×
            nn.Linear(3*d_model, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, num_labels)
        )
        self.cross_attn = nn.MultiheadAttention(
            embed_dim=d_model,  # 输入/输出维度
            num_heads=nhead,  # 与Transformer头数一致
            dropout=dropout,  # 复用dropout率
            batch_first=True  # 批处理优先格式
        )

    def forward(self, vec, labels=None):
        batch_size, device = vec.size(0), vec.device

        # 拆分输入向量
        e_src = vec[:, :32]  # 源基因表达特征
        g_src = vec[:, 32:32 + 13805]  # 源GO特征
        e_tgt = vec[:, 32 + 13805:32 + 13805 + 32]  # 目标基因表达特征
        g_tgt = vec[:, 32 + 13805 + 32:]  # 目标GO特征

        # 处理源基因表达→32个独立token
        src_e = self.expr_proj(e_src).view(batch_size, 32, -1)  # [B, 32, d_model]
        # 处理源GO→64个token
        src_g = self.go_down(g_src).view(batch_size, 64, -1)  # [B, 64, d_model]
        # 拼接源特征
        src = torch.cat([src_e, src_g], dim=1)  # [B, 96, d_model]

        # 处理目标基因表达→32个独立token
        tgt_e = self.expr_proj(e_tgt).view(batch_size, 32, -1)  # [B, 32, d_model]
        # 处理目标GO→64个token
        tgt_g = self.go_down(g_tgt).view(batch_size, 64, -1)  # [B, 64, d_model]
        # 拼接目标特征
        tgt = torch.cat([tgt_e, tgt_g], dim=1)  # [B, 96, d_model]

        # 合并源+目标序列
        seq = torch.cat([src, tgt], dim=1)  # [B, 192, d_model]
        # 添加位置编码
        seq = seq + self.pos_embed.to(device)

        # Transformer编码
        encoded = self.encoder(seq)  # [B, 192, d_model]
        src_seq = encoded[:, :96, :]  # 源基因序列 [B, 96, d_model]
        tgt_seq = encoded[:, 96:, :]  # 目标基因序列 [B, 96, d_model]

        # 融合原有特征（拼接后降维）
        src_embed = src_seq.mean(dim=1)  # 源基因全局特征 [B, d_model]
        tgt_embed = tgt_seq.mean(dim=1)  # 目标基因全局特征 [B, d_model]


        # 交叉注意力：目标基因作为query，源基因作为key和value
        tgt_query_src, _ = self.cross_attn(
            query=tgt_seq,
            key=src_seq,
            value=src_seq
        )  # [B, 96, d_model]
        src_query_tgt, _ = self.cross_attn(
            query=src_seq,
            key=tgt_seq,
            value=tgt_seq
        )  # [B, 96, d_model]
        tgt_to_src_embed = tgt_query_src.mean(dim=1)  # [B, d_model]
        src_to_tgt_embed = src_query_tgt.mean(dim=1)
        diff = (tgt_embed - src_embed)
        # 拼接三个特征向量
        # fused = torch.cat([src_embed, tgt_embed, diff], dim=-1)  # [B, 3*d_model]
        # fused = torch.cat([diff], dim=-1)  # [B, 1*d_model]
        # fused = torch.cat([tgt_to_src_embed,src_to_tgt_embed], dim=-1)  # [B, 2*d_model]
        # fused = torch.cat([tgt_to_src_embed, src_to_tgt_embed, src_embed, tgt_embed, diff], dim=-1)  # [B, 2*d_model]
        # fused = torch.cat([src_embed, tgt_embed], dim=-1)  # [B, 2*d_model]
        fused = torch.cat([tgt_to_src_embed, src_to_tgt_embed, diff], dim=-1)  # [B, 3*d_model]



        # 调整分类器输入维度
        logits = self.classifier(fused)  # 需同步修改classifier第一层输入尺寸
        # 计算损失
        loss = None
        if labels is not None:
            loss = torch.nn.CrossEntropyLoss()(logits, labels)

        return {"loss": loss, "logits": logits}


# print(model)
# dummy_input = torch.randn(1, 27674)
# traced_model = torch.jit.trace(model, dummy_input)
# traced_model.save("EGO_former.pt")
# 修改 compute_metrics 函数
def compute_metrics(eval_pred):
    logits, labels = eval_pred
    print(eval_pred)
    # 处理可能的元组嵌套
    if isinstance(logits, tuple):
        # print(logits)
        logits = logits[0]

    # 计算预测结果
    preds = np.argmax(logits, axis=-1)


    # 计算概率（二分类）
    if logits.shape[1] == 2:
        probs = softmax(logits, axis=-1)[:, 1]
    else:
        probs = softmax(logits, axis=-1)
    acc = accuracy_score(labels, preds)
    precision_macro = precision_score(labels, preds, average='macro')
    neg_recall = recall_score(labels, preds, pos_label=0)
    pos_recall = recall_score(labels, preds, pos_label=1)
    g_mean = (neg_recall * pos_recall) ** (1 / 2)
    f1 = f1_score(labels, preds, average='macro')
    # cm = confusion_matrix(labels, preds)
    auroc = roc_auc_score(labels, probs)
    aupr = average_precision_score(labels, probs)
    mcc = matthews_corrcoef(labels, preds)

    return {
        "A_accuracy": acc,
        "B_precision":precision_macro,
        "C-G-mean":g_mean,
        "D_f1": f1,
        "E_mcc": mcc,
        "F_auroc": auroc,
        "G_aupr":aupr,
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--txt_input", required=True, help="nature_amps.csv")
    parser.add_argument("--array_input", required=True, help="Artifical_amps.csv")
    parser.add_argument("--output_dir", required=True, help="Artifical_amps.csv")
    parser.add_argument("--count_fold", required=True, help="Artifical_amps.csv")
        # parser.add_argument("--output_dir", required=True, help="Artifical_amps.csv")
    args = parser.parse_args()


    # output_dir = "./conclude_new_three"
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    expr_go = np.load(args.array_input).astype(np.float32)
    X, y = [], []
    val_idx = []                          # 记录验证集样本下标
    # with open('../../hsa_bar_image_5k.txt') as f:
    with open(args.txt_input) as f:
        for idx, line in enumerate(f):
            parts = line.strip().split('\t')
            src_idx   = int(parts[2])
            tgt_idx   = int(parts[3])
            label     = 1 if parts[4] == '"regulation"' else 0
            sample_tp = parts[5]

            vec = np.concatenate([expr_go[src_idx], expr_go[tgt_idx]])
            X.append(vec.astype(np.float32))
            y.append(label)

            # 划分验证集
            if sample_tp == args.count_fold:
                val_idx.append(idx)

    X = np.stack(X)

    y = np.array(y, dtype=np.int64)

    # ---------- 3. 训练/验证划分 ----------
    val_idx = np.array(val_idx)
    train_mask = np.ones(len(y), dtype=bool)
    train_mask[val_idx] = False
    X_train, y_train = X[train_mask], y[train_mask]
    X_val,   y_val   = X[val_idx],   y[val_idx]
    print(f"train={len(y_train)}, val={len(y_val)}, feat_dim={X.shape[1]}")

    # ---------- 4. Dataset ----------
    
    train_ds = VecDataset(X_train, y_train)
    val_ds   = VecDataset(X_val,   y_val)

    model = EGO_former(d_model=128, nhead=8, num_layers=4)

    args = TrainingArguments(
    output_dir=output_dir,
    num_train_epochs=30,
    per_device_train_batch_size=64,
    per_device_eval_batch_size=128,
    learning_rate=1e-6,
    weight_decay=0.01,
    lr_scheduler_type="cosine",
    warmup_ratio=0.1,
    # max_grad_norm=1.0 ,       # 正确写法
    fp16=torch.cuda.is_available(),
    evaluation_strategy="epoch",
    save_strategy="epoch",
    save_total_limit=2,
    load_best_model_at_end=True,
    metric_for_best_model="A_accuracy",
    report_to="none",
    # dataloader_num_workers=4,
    seed=42,
    )
    trainer = Trainer(
        model=model,
        args=args,
        train_dataset=train_ds,
        eval_dataset=val_ds,
        compute_metrics=compute_metrics,
        tokenizer=None,
    )
    trainer.train()

    best = os.path.join(output_dir, "best_seed42")
    trainer.save_model(best)
    print("训练完成，已保存到", best)