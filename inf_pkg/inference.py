

from pickle import load
from utils.datasets import ProteinDataset
from argparse import ArgumentParser
from transformers import (
    EsmTokenizer,
    EsmForSequenceClassification,
    EsmConfig,
)
from peft import (
    LoraConfig,
    PeftModel,
    get_peft_model,
)
from tqdm import tqdm
from torch.utils.data import DataLoader
from torch import stack, argmax
from numpy import var, concatenate
from pandas import DataFrame
import torch

print(torch.cuda.is_available())
DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
MODEL = 'esm2_8m_base'
peft_model_id = 'cycformer_lora_80_qk/checkpoint-157636'

torch.manual_seed(42)

def main():
    print('DEV: ', DEVICE)
    args = parse_args()
    tokenizer = EsmTokenizer.from_pretrained('facebook/esm2_t6_8M_UR50D')
    with open('ko_model_id_map.pickle', 'rb') as f:
        mapper = load(f)

    dataset = args.proteins
    print(dataset)
    if 'fasta_file' in args:
        val_proteins = ProteinDataset.from_fastx(args.fasta_file, 
                                                tokenizer=tokenizer, 
                                                mapper=mapper,
                                                )
    else:
        val_proteins = ProteinDataset(dataframe=None,
                                      return_seqs=True,
                                      csv_file=None, 
                                      tokenizer=tokenizer, 
                                      mapper=mapper,
                                      nsamples=-1,
                                      seq_column='sequence',
                                      label_column='label')
    
    model = EsmForSequenceClassification.from_pretrained(MODEL,
                                                        )
    
    lora_model = PeftModel.from_pretrained(model, peft_model_id)
    lora_model.merge_and_unload()
    model = lora_model
    model.to(DEVICE)

    # DROPOUT FUNCTIONALITY. DONT TURN ON...YET
    #monte_carlo_dropout(model, val_proteins, 'protein_annotations')

    predict(model, val_proteins, args.annotations)

def predict(model, dataset, save_filename):
    loader = DataLoader(dataset, shuffle=False, batch_size=1)
    predictions = []
    labels = []
    seqs = []
    embeddings = []
    for sample in tqdm(loader):
        outputs = model(sample['input_ids'].to(DEVICE).squeeze(0), output_hidden_states=True)

        embedding = outputs.hidden_states[-1][:, 0, :].detach().cpu().numpy()
        embeddings.append(embedding)

        y_hat = argmax(outputs.logits, dim=-1)
        y_hat = dataset.id2label[y_hat.detach().cpu().numpy()[0]]
        
        if 'label' in sample:
            labels.append(dataset.id2label[sample['label'].numpy()[0]])
        predictions.append(y_hat)
        seqs.append(sample['sequence'][0])

    embeddings = concatenate(tuple(embeddings), axis=0)
    annotations = {}
    annotations['sequence'] = seqs
    annotations['predictions'] = predictions
    df = DataFrame(annotations)
    df.to_csv(save_filename+'.csv')
    print('done')
    return df

def monte_carlo_dropout(model, dataset, filename, rounds=50):
    model.train()
    loader = DataLoader(dataset, shuffle=False, batch_size=1)
    seq_annotations = {}
    annot = []
    seqs = []
    v = []
    labels = []
    model = model.to(DEVICE)
    for sample in tqdm(loader):
        preds = stack([argmax(model(sample['input_ids'].to(DEVICE).squeeze(0)).logits, dim=-1) 
                       for _ in range(rounds)]).squeeze(-1).cpu().numpy()
        variance = var(preds)
        preds = [dataset.id2label[idx] 
                 for idx in preds]
        if 'label' in sample:
            labels.append(dataset.id2label[sample['label'].numpy()[0]])
        seqs.append(sample['sequence'][0])
        annot.append(preds)
        v.append(variance)

    seq_annotations['sequence'] = seqs
    seq_annotations['ANNOT'] = annot
    seq_annotations['variance'] = v
    if 'label' in sample:
        seq_annotations['label'] = labels
    df = DataFrame(seq_annotations)
    df.to_csv(filename+'.csv')
    print(df)
    return df

def parse_args():
    # parse arguments
    parser = ArgumentParser()
    parser.add_argument('--fasta_file')
    parser.add_argument('--annotations')
    parser.add_argument('--proteins')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()