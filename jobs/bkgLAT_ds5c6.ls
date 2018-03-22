don't submit this one, it has missing files from the batchSplit fkup

./wave-skim -n -r 5 112 -p /global/projecta/projectdirs/majorana/users/wisecg/bkg/skim .
./wave-skim -n -r 5 113 -p /global/projecta/projectdirs/majorana/users/wisecg/bkg/skim .

./lat.py -b -c -r 5 120 -p /global/projecta/projectdirs/majorana/users/wisecg/bkg/split/splitSkimDS5_120.root ./latSkimDS5_120_0.root
./lat.py -b -c -r 5 106 -p /global/projecta/projectdirs/majorana/users/wisecg/bkg/split/splitSkimDS5_106_2.root ./latSkimDS5_106_2.root