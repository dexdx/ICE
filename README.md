# What's ICE?

ICE (Indirect Cointegration Estimation) is a new, (very) simple and indirect method for estimating cointegrating relationships in multivariate stochastic processes. The method exploits Johansen's benchmark representation of a cointegrated system [(reduced VECM)](https://www.sciencedirect.com/science/article/pii/0165188988900413) to estimate the cointegrating matrix via a VAR estimation. Because of that, it is a *generic* approach in which any VAR estimator can be applied, and as a consequence, is also applicable in high dimensions and/or sparse settings by accordingly choosing a regularized estimator.

<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

