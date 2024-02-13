import wandb
from sklearn.model_selection import train_test_split



wandb.login()

run = wandb.init(project='cycle_pred')

artifact = wandb.Artifact(name='cycle-pred',
                          type='dataset')

