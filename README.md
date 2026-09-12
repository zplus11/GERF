<a href="https://doi.org/10.5281/zenodo.20706333"><img src="https://zenodo.org/badge/1197217781.svg" alt="DOI"></a>

<img align="left" src="img/gerf.png" width="200"> is a Mathematica paclet that solves nonlinear integral order partial differential equations using the GERF (generalized exponential rational function) expansion technique. For more info on the method, see [10.1140/epjp/i2018-11984-1](https://doi.org/10.1140/epjp/i2018-11984-1). To learn more about the paclet, visit the homepage at [Wolfram repository](https://resources.wolframcloud.com/PacletRepository/resources/Taggar/GERF/).

# Installation

is as easy as running

```mathematica
PacletInstall["Taggar/GERF"]
```

in your local or cloud notebook. Then, load it with

```mathematica
<<Taggar`GERF`
```

If you decide to use the paclet for research, kindly consider citing it as follows:

Taggar, N. (2026). GERF: A Mathematica paclet for the generalized exponential rational function method (1.2.1). Zenodo. [https://doi.org/10.5281/zenodo.20706333](https://doi.org/10.5281/zenodo.20706333).

# Usage

Let

```mathematica
burgers = D[u[x, t], t] + u[x, t] * D[u[x, t], x] - \[Nu] * D[u[x, t], x, x] == 0
```

be the given equation (this is the Burgers' equation in (1+1)-dimensions). Then, use GERFSolve as follows:

```mathematica
sol = GERFSolve[burgers, u[x, t]]
```

which returns the following output:

![Solutions](img/sols.png)

Pick any of these and plot it for appropriate values:

```mathematica
Plot3D[
    u[x, t] /. sol[[3]] /. {(* parameters *)},
    {x, -4, 4}, {t, 0, 4},
    PlotRange -> All]
```
![Graph](img/plot.png)

Furthermore, just the balance constant of an equation can be computed. For example:

```mathematica
BalanceConstant[burgers, u[x, t]]
(* 1 *)
```

A full demonstration is available on the Wolfram repository page of this paclet, see [GERF](https://resources.wolframcloud.com/PacletRepository/resources/Taggar/GERF).

# Version log

**Version 1.3.0,** *05 September, 2026* &mdash; some code refactoring and introducing `BalanceConstant` in the public API.

**Version 1.2.1,** *24 June, 2026* &mdash; additional examples added in documentation, and citation information also added.

**Version 1.2.0,** *15 June, 2026* &mdash; support for systems of equations.

**Version 1.1.0,** *on 30 May, 2026* &mdash; support for fractional ordered equations.

**Version 1.0.0,** *on 31 March, 2026* &mdash; initial upload.

# Contributions

are welcome.
