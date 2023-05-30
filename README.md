# binSIDH and terSIDH

A proof-of-concept implementation of binSIDH, terSIDH, and their hybrid variants, as proposed in New SIDH Countermeasures for a More Efficient Key Exchange by Andrea Basso and Tako Boris Fouotsa.

# How to use

The implementation of binSIDH and terSIDH is in `bin-terSIDH.sage`, while that of the hybrid variants binSIDH<sup>hyb</sup> and terSIDH<sup>hyb</sup> is in `bin-terSIDH--hybrid.sage`. The remaining files are sourced from the [Kummer Isogeny library](https://github.com/jack4818/KummerIsogeny) by Giacomo Pope and the [FESTA SageMath implementation](https://github.com/FESTA-PKE/FESTA-SageMath), on which these implementations are based on.

The four protocols can be run with Sage, with the following arguments:

```
sage bin-terSIDH.sage [--128, --192, --256] [--bin, --ter]
```

By default, the 128-bit security parameters are selected. To access other security levels:
  - The flag `--192` selects the parameters aiming for 192-bit security 
  - The flag `--256` selects the parameters aiming for 256-bit security 

Similarly, by default, the ternary version is run. The binary variant is selected with the `--bin` flag.
