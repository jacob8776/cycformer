python train.py --model_save_dir=cycformer --run=cycformer --project=protein-cycle-pred --artifact=cycle-pred/protein-cycle-pred/protein-sequences:v0 --trainseqs=train_seqs --valseqs=val_seqs --model=cycformer/checkpoint-76146 --tokenizer=facebook/esm2_t6_8M_UR50D --label=CYCLE --inputs=SEQS 