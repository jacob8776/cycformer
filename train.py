
import numpy as np
import torch
import wandb
import argparse
from evaluate import load
from transformers import (
    Trainer,
    TrainingArguments,
    EsmTokenizer,
    EsmForSequenceClassification, 
)
from utils.datasets import ProteinDataset

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

    print(f'PROJECT:         {PROJECT_NAME}')
    print(f'TRAINING DATA:   {TRAIN_DATA_FILENAME}')
    print(f'VALIDATION DATA: {VAL_DATA_FILENAME}')
    print(f'ARTIFACT:        {ARTIFACT_FILENAME}')
    print(f'MODEL:           {MODEL}')
    print(f'TOKENIZER:       {TOKENIZER_NAME}')


    # download data and get filenames
    run = wandb.init(project=PROJECT_NAME, mode='disabled')
    artifact = run.use_artifact(ARTIFACT_FILENAME, type='dataset')
    artifact.download()

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
                      tokenizer=tokenizer,
                      compute_metrics=compute_metrics,
                    )

    trainer.train()
    trainer.evaluate()
    wandb.finish()



def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)

    precision = precision_metric.compute(predictions=predictions, references=labels, average='weighted')
    recall = recall_metric.compute(predictions=predictions, references=labels, average='weighted')
    acc = accuracy.compute(predictions=predictions, references=labels)
    mcc = matthews_metric.compute(references=labels, predictions=predictions)    
    f1 = f1_metric.compute(references=labels, predictions=predictions, average='weighted')
    return {"accuracy": acc, "MCC": mcc, "F1": f1, 'precision': precision, 'recall': recall}

def parse_args():
    # parse arguments
    parser = argparse.ArgumentParser()
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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
    
                        

