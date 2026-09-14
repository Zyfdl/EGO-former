# EGO-former

Predicting the interactions between transcription factors and their target genes is not only crucial for understanding the regulation of gene expression in organisms ,but also helpful for the diagnosis and treatment of diseases.To infer new transcriptional regulatory relationships,we built a new Transformer-based model EGO-former for transcription regulation prediction,with regarding Single Cell Expression Atlas
and GO terms as the sample data.EGO-former is a Encoder-only model that combined a bidirectional cross attention module to extract the interaction features between transcription factors and target genes with a difference of the target gene to the transcription factor encoding pooled vectors to capture more feature informations. Following 5-fold cross experiments,our model not only achieved a prediction accuracy of 0.8771 for transcription factors and their target genes that was outperforming classical models such as
Convolutional Neural Networks and Graph Neural Networks,but also demonstrated superior performance across other evaluation metrics.

![GO_forme](.\EGO_former.png)

## Running

Dependencies (with python >= 3.9).Example commands to install the dependencies in a new conda environment.

```
conda create --name tftarget python=3.9
conda activate tftarget
pip install -r requirements.txt --extra-index-url https://download.pytorch.org/whl/cu128
nohup bash run_EGO_former.sh &
```

