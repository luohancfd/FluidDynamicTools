## 0DRate calculation

1. Set the following flags
    - IREAC = 2
    - IRELAX = 0
    - IMFS = 1 if you use MF-DSMC model
2. Create `DS1VD.in` and `input.txt` like the ones in folder `Example/0DRate`. Remember to change `IRM` to match the reaction you want as well as the appropriate mole fraction of the species

## 0DRelax calculation
1. Set the following flags
    - IREAC = 0
    - IRELAX = 0
2. Create `DS1VD.in` and `input.txt` like the ones in folder `Example/0DRelax`. `IRM` has to be `0` to turn off all reactions
3. Track `RELAX.TXT`, it should contain information needed