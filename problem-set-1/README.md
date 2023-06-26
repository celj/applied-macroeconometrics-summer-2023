# Problem Set 1

File structure:

```sh
.
├── README.md             <- this file
├── chord.m               <- chord method function
├── data                  <- data directory
│   └── LM_grossflows.xls <- data file
├── gauss_seidel.m        <- gauss seidel method function
├── jacobi.m              <- jacobi method function
├── labor_ex.m            <- labor exercise
├── linear_ex.m           <- linear exercise
├── newton.m              <- newton method function
├── nonlinear_ex.m        <- nonlinear exercise
├── setup.m               <- settings file
└── src.m                 <- main file
```

Only `src.m` is required to run the code. It calls `setup.m` and every `*_ex.m` file.

> The other files are helper functions or data files.
