# Resolved top tagger training using xgboost

Install the prerequisites:

- `pip install weaver-core xgboost`

To run the training:

- update the file path `train_val_files` in the script
  - comment out the line `train_val_files = train_val_files[:2]` to use all files for the final training, or increase the number of files
- run the script:
  - `python xgb_train_resTop --train`
