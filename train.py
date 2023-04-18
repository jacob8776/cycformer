
import numpy as np
import torch
import pandas as pd
from tqdm import tqdm
from transformers import AutoTokenizer, Trainer, TrainingArguments
from transformers.modeling_outputs import SequenceClassifierOutput
from transformers import AutoModelForSequenceClassification
from transformers import EsmModel, EsmTokenizer, EsmForSequenceClassification
from transformers.trainer_utils import EvalLoopOutput
from torch.utils.data import DataLoader
from torch.optim import AdamW
from transformers import get_scheduler
from datasets import load_metric
import evaluate
import wandb
import argparse
from evaluate import load, save
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import numpy as np
import pandas as pd
from Bio import SeqIO


torch.manual_seed(42)

accuracy = load("accuracy")
matthews_metric = load("matthews_correlation")
f1_metric = load("f1")
precision_metric = load("precision")
recall_metric = load('recall')

def main():
    # main logic
    args = parse_args()
    output_dir = args.model_save_dir
    run_name = args.run
    PROJECT_NAME = args.project
    TRAIN_DATA_FILENAME = args.trainseqs
    VAL_DATA_FILENAME = args.valseqs
    ARTIFACT_FILENAME = args.artifact
    MODEL = args.model
    TOKENIZER_NAME = args.tokenizer
    LABEL_NAME = args.label
    INPUT_NAME = args.inputs
    ONLY_EVAL = args.eval
    FASTA_FILE = args.fasta_file

    print(f'PROJECT:         {PROJECT_NAME}')
    if ONLY_EVAL:
        print('INFERENCE MODE')
        print(FASTA_FILE)
    else:
        print(f'TRAINING DATA:   {TRAIN_DATA_FILENAME}')
        print(f'VALIDATION DATA: {VAL_DATA_FILENAME}')
        print(f'ARTIFACT:        {ARTIFACT_FILENAME}')
        print(f'MODEL:           {MODEL}')
        print(f'TOKENIZER:       {TOKENIZER_NAME}')


    # download data and get filenames
    run = wandb.init(project=PROJECT_NAME)
    artifact = run.use_artifact(ARTIFACT_FILENAME, type='dataset')
    artifact_dir = artifact.download()

    data_folder = ARTIFACT_FILENAME.split('/')[-1]
    TRAIN_DATA_FILENAME = 'artifacts'+ '/' + data_folder + '/' + TRAIN_DATA_FILENAME
    VAL_DATA_FILENAME = 'artifacts' + '/' + data_folder + '/' + VAL_DATA_FILENAME

    # set up tokenizer
    tokenizer = EsmTokenizer.from_pretrained(TOKENIZER_NAME)

    # get datasets set up
    train_proteins = ProteinDataset(csv_file=TRAIN_DATA_FILENAME, 
                                    tokenizer=tokenizer,
                                    label_column=LABEL_NAME,
                                    seq_column=INPUT_NAME,
                                    dataframe=None,
                                    )
    val_proteins = ProteinDataset(csv_file=VAL_DATA_FILENAME, 
                                tokenizer=tokenizer,
                                label_column=LABEL_NAME,
                                seq_column=INPUT_NAME,
                                nsamples=3,
                                mapper=train_proteins.label_dict,
                                dataframe=None,
                                )
    
    # check if we are just doing predictions
    if args.eval:
        mapper = train_proteins.label_dict
        val_proteins = ProteinDataset.from_fastx(FASTA_FILE, tokenizer=tokenizer, mapper=mapper)

    # init model
    model = EsmForSequenceClassification.from_pretrained(MODEL, 
                                                         num_labels=train_proteins.nu_labels,
                                                         hidden_dropout_prob=0.1,
                                                         attention_probs_dropout_prob=0.1,
                                                        )

    # model training
    model.train()
    print(f'# of TRAIN SEQ:   {len(train_proteins)}')
    print(f'# of VAL SEQ:     {len(val_proteins)}')
    print(f'# OF LABELS:      {train_proteins.nu_labels}')

    tr_args = TrainingArguments(output_dir,
                                evaluation_strategy='steps',
                                save_strategy='steps',
                                learning_rate=3e-4,
                                num_train_epochs=3,
                                #max_steps=500000,
                                logging_steps=1,
                                auto_find_batch_size=True,
                                report_to='wandb',
                                run_name=run_name,
                                save_steps=30000,
                                eval_steps=10000,
                                eval_accumulation_steps=100,
                            )

    trainer = Trainer(model,
                      tr_args,
                      train_dataset=train_proteins,
                      eval_dataset=val_proteins,
                      #tokenizer=tokenizer,
                      compute_metrics=compute_metrics,
                    )

    # eval logic
    if args.eval:
        monte_carlo_dropout(model, val_proteins, filename='cycle_annotations', rounds=50)
    else:
        trainer.train()
        trainer.evaluate()
    wandb.finish()



"""
For loading dataframes containing protein sequences with a target ID
This serves as a general dataset for protein sequence classification

Arguments:
csv_file: the CSV file name that will be used as a pandas DataFrame
tokenizer: the HuggingFace tokenizer for processing sequences
label_column: the column name in the DataFrame corresponding to the targets
seq_column: the column name in the DataFrame corresponding to the proteins
"""
class ProteinDataset(Dataset):
    def __init__(self, dataframe, csv_file, tokenizer, label_column, seq_column, nsamples=-1, mapper=None):
        self.tokenizer = tokenizer
        data = dataframe
        if dataframe is None:
            print('loading file: ', csv_file)
            data = pd.read_csv(csv_file)
        
        self.dataframe = data
        if nsamples > 0:
            data = data.groupby(label_column).sample(n=nsamples).reset_index()
        print('done.')
        self.label_dict = mapper
        if label_column is not None:
            self.labels = data[label_column]
            self.nu_labels = len(np.unique(self.labels))
            if mapper == None:
                label_dict = {}
                label_list = np.unique(self.labels)
                for i in range(len(label_list)):
                    label_dict[label_list[i]] = i
                self.label_dict = label_dict
            #else: self.label_dict = mapper
            self.label_name = label_column

        self.id2label = dict([(k, v) for v, k in self.label_dict.items()])
        self.seqs = data[seq_column]

    """
    This should work for a correctly formatted fasta file!
    """
    @classmethod
    def from_fastx(cls, fasta_filename, tokenizer, seq_column='SEQS', label_column=None, nsamples=-1, mapper=None):
        # read fasta file into dataframe
        df = None
        file_suffix = fasta_filename[-5:]
        print('File type: ', file_suffix)
        with open(fasta_filename) as fasta_file:  # Will close handle cleanly
            identifiers = []
            seqs = []
            for seq_record in SeqIO.parse(fasta_file, file_suffix):  # (generator)
                identifiers.append(seq_record.id)
                seqs.append(str(seq_record.seq))
            df = pd.DataFrame({'IDs': identifiers,'SEQS': seqs})

            
            df = df.sample(n=100).reset_index() # test

            return cls(dataframe=df, 
                    csv_file=None, 
                    tokenizer=tokenizer, 
                    seq_column=seq_column, 
                    label_column=label_column,
                    mapper=mapper
                    )
            
    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, idx):
        sample = self.tokenizer(self.seqs[idx], #.replace(' ', ''), 
                                return_tensors='pt', 
                                padding=True,
                                truncation=True,
                                )
        #sample['input_ids'] = sample['input_ids'].squeeze(0)
        #print(sample)
        if hasattr(self, 'labels'):
            sample['label'] = self.label_dict[self.labels[idx]]
        sample['sequence'] = self.seqs[idx]
        return sample

def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)

    precision = precision_metric.compute(predictions=predictions, references=labels, average='weighted')
    recall = recall_metric.compute(predictions=predictions, references=labels, average='weighted')
    acc = accuracy.compute(predictions=predictions, references=labels)
    #probs = torch.nn.functional.softmax(torch.Tensor(logits), dim=1)
    mcc = matthews_metric.compute(references=labels, predictions=predictions)    
    f1 = f1_metric.compute(references=labels, predictions=predictions, average='weighted')
    return {"accuracy": acc, "MCC": mcc, "F1": f1, 'precision': precision, 'recall': recall}

def parse_args():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file')
    parser.add_argument('--model_save_dir')
    parser.add_argument('--run')
    parser.add_argument('--project')
    parser.add_argument('--artifact')
    parser.add_argument('--trainseqs')
    parser.add_argument('--valseqs')
    parser.add_argument('--model')
    parser.add_argument('--tokenizer')
    parser.add_argument('--label')
    parser.add_argument('--inputs')
    parser.add_argument('--eval', action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    return args

# TODO
def get_datasets(args):
    pass

# TODO
def download_data(args):
    pass

# TODO
# Add computations of variance per class over all sequences
# can we also get attention map visuals?
def monte_carlo_dropout(model, dataset, filename, rounds=50):
    model.train()
    loader = DataLoader(dataset, shuffle=False, batch_size=1)
    seq_annotations = {}
    annot = []
    seqs = []
    v = []

    for sample in tqdm(loader):

        preds = torch.stack([torch.argmax(model(sample['input_ids'].squeeze(0)).logits, dim=-1) 
                             for _ in range(rounds)]).squeeze(-1).numpy()
        
        variance = np.var(preds)

        preds = [dataset.id2label[idx] 
                 for idx in preds]
        
        seqs.append(sample['sequence'][0])
        annot.append(preds)
        v.append(variance)

    seq_annotations['sequence'] = seqs
    seq_annotations['ANNOT'] = annot
    seq_annotations['variance'] = v
    df = pd.DataFrame(seq_annotations)
    df.to_csv(filename+'.csv')
    print(df)
    return df
        
        
        

if __name__ == '__main__':
    main()
    
                        

