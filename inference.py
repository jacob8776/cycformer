# TODO

# add logic for running inference on protein sequences here

# need ProteinDataset class

# need to save KO IDs and load them in for inference step

# 

from pickle import load
from utils.datasets import ProteinDataset
from argparse import ArgumentParser
from transformers import (
    EsmTokenizer,
    EsmForSequenceClassification, 
)
from tqdm import tqdm
from torch.utils.data import DataLoader
from torch import stack, argmax
from numpy import var
from pandas import DataFrame

def main():
    
    args = parse_args()
    tokenizer = EsmTokenizer.from_pretrained('facebook/esm2_t6_8M_UR50D')
    with open('ko_model_id_map.pickle', 'rb') as f:
        mapper = load(f)

    val_proteins = ProteinDataset.from_fastx(args.fasta_file, 
                                             tokenizer=tokenizer, 
                                             mapper=mapper,
                                             )
    print(val_proteins.id2label)
    model = EsmForSequenceClassification.from_pretrained('cycformer/checkpoint-76146',
                                                          hidden_dropout_prob=0.1,
                                                          attention_probs_dropout_prob=0.1,
                                                        )
    monte_carlo_dropout(model, val_proteins, 'protein_annotations')


def monte_carlo_dropout(model, dataset, filename, rounds=50):
    model.train()
    loader = DataLoader(dataset, shuffle=False, batch_size=1)
    seq_annotations = {}
    annot = []
    seqs = []
    v = []

    for sample in tqdm(loader):
        preds = stack([argmax(model(sample['input_ids'].squeeze(0)).logits, dim=-1) 
                       for _ in range(rounds)]).squeeze(-1).numpy()
        variance = var(preds)
        preds = [dataset.id2label[idx] 
                 for idx in preds]
        seqs.append(sample['sequence'][0])
        annot.append(preds)
        v.append(variance)

    seq_annotations['sequence'] = seqs
    seq_annotations['ANNOT'] = annot
    seq_annotations['variance'] = v
    df = DataFrame(seq_annotations)
    df.to_csv(filename+'.csv')
    print(df)
    return df

def parse_args():
    # parse arguments
    parser = ArgumentParser()
    parser.add_argument('--fasta_file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()