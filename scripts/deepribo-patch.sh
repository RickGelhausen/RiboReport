#!/usr/bin/env bash

trap ctrl_c INT

function ctrl_c() {
echo "** Trapped CTRL-C"
exit;
}

sed -i 's/model.load_state_dict(torch.load(model_name, map_location=device))/model.load_state_dict(torch.load(model_name, map_location=device),strict=False)/g' ../../tools/DeepRibo/src/DeepRibo.py
sed -i 's/sys.exit(ParseArgs())/ParseArgs()/g' ../../tools/DeepRibo/src/DeepRibo.py
exit;
