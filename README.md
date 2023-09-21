# Celestial singularities

We present Maple code for the computations of singularities of celestial surfaces
that are represented in terms of their parametric types.

For the definitions and details we refer to the article

[Self-intersections of surfaces that contain two circles through each point](https://arxiv.org/abs/?).

See also the documentation in [ParametricTypes.m](https://github.com/niels-lubbes/celestial-singularities/blob/master/ParametricTypes.m).

A parametric type corresponds to a birational map from P1xP1
to a surface X that lies in S3, where P1xP1 denote the
product of the projective line with itself and S3 denotes the projectivized unit 3-sphere.
After taking an affine chart, we obtain a map from complex plane C2 to X.
We compute the components of the preimage in C2 of the singular components in X
such that a general complex point of this component in X
has at least two preimages in the parameter domain C2.

Download [Maple](https://www.maplesoft.com/products/maple/)
and
copy paste the code at
[ParametricTypes.m](https://github.com/niels-lubbes/celestial-singularities/blob/master/ParametricTypes.m)
into a Maple session.
The code has been tested in Maple 13.
To run from a Linux shell we can use the following command:

```
maple < ParametricTypes.m
```

Output

```
parametric type="D3a" value=[1, -1, 0, 0, 1, 0, 3/2, 0]
        bideg=[0, 1], sing=[1], mult=1, factor=t+1/29*RootOf(_Z^2+145)
        bideg=[0, 1], sing=[1], mult=1, factor=t-1/29*RootOf(_Z^2+145)
        bideg=[1, 0], sing=[1], mult=1, factor=s-RootOf(3*_Z^2-4*_Z+3)
        bideg=[1, 0], sing=[1], mult=1, factor=-4/3+RootOf(3*_Z^2-4*_Z+3)+s
        bideg=[2, 2], sing=[1], mult=1, factor=s^2*t^2-72/235*t*s^2-29/235*s^2-538/235*s*t^2+86/235*s+223/235*t^2-72/235*t-41/235
parametric type="D3b" value=[1/10, 1/2, 0, 0, 1/2, 0, 1/2, 0]
        bideg=[0, 1], sing=[1], mult=1, factor=t+1/2*RootOf(_Z^2+2)
        bideg=[0, 1], sing=[1], mult=1, factor=t-1/2*RootOf(_Z^2+2)
        bideg=[2, 2], sing=[1], mult=1, factor=s^2*t^2+1525/837*t*s^2+20/93*s^2-10/27*s*t^2+110/93*t*s+250/837*s+47/27*t^2+1525/837*t-130/837
        bideg=[1, 0], sing=[1], mult=1, factor=s-RootOf(63*_Z^2+10*_Z+63)
        bideg=[1, 0], sing=[1], mult=1, factor=10/63+RootOf(63*_Z^2+10*_Z+63)+s
parametric type="D3c" value=[3/10, 1/2, 0, 0, 1/2, 0, 1/2, 0]
        bideg=[1, 0], sing=[1], mult=1, factor=s-RootOf(67*_Z^2+30*_Z+67)
        bideg=[1, 0], sing=[1], mult=1, factor=s+30/67+RootOf(67*_Z^2+30*_Z+67)
        bideg=[2, 2], sing=[1], mult=1, factor=s^2*t^2+1225/87*t*s^2+220/29*s^2-10*s*t^2+910/29*t*s+250/29*s+21*t^2+1225/87*t-70/29
        bideg=[0, 1], sing=[1], mult=1, factor=t+1/2*RootOf(_Z^2+2)
        bideg=[0, 1], sing=[1], mult=1, factor=t-1/2*RootOf(_Z^2+2)
parametric type="D3d" value=[1, -4, 0, 0, 1, 0, 0, -1]
        bideg=[2, 2], sing=[1], mult=1, factor=s^2*t^2-240/331*t*s^2+171/331*s^2-216/331*s-216/331*s*t^2+395/331*t^2+16/331*t+107/331
        bideg=[0, 1], sing=[1], mult=1, factor=t-RootOf(1+_Z^2)
        bideg=[0, 1], sing=[1], mult=1, factor=RootOf(1+_Z^2)+t
        bideg=[1, 0], sing=[1], mult=1, factor=s-RootOf(9*_Z^2-8*_Z+9)
        bideg=[1, 0], sing=[1], mult=1, factor=s-8/9+RootOf(9*_Z^2-8*_Z+9)
parametric type="D3e" value=[1/5, 1/2, 0, 0, 1/2, 0, 1/2, 0]
        bideg=[2, 2], sing=[1], mult=1, factor=s^2*t^2+11300/3751*t*s^2+3220/3751*s^2-40/31*s*t^2+15360/3751*t*s+4000/3751*s+111/31*t^2+11300/3751*t-1620/3751
        bideg=[1, 0], sing=[1], mult=1, factor=40/129+RootOf(129*_Z^2+40*_Z+129)+s
        bideg=[1, 0], sing=[1], mult=1, factor=s-RootOf(129*_Z^2+40*_Z+129)
        bideg=[0, 1], sing=[1], mult=1, factor=t+1/2*RootOf(_Z^2+2)
        bideg=[0, 1], sing=[1], mult=1, factor=t-1/2*RootOf(_Z^2+2)
parametric type="D4a" value=[1, 1, 0, 0, 1, 0, 1, 0]
        bideg=[0, 1], sing=[1], mult=1, factor=t-1/5*RootOf(_Z^2+5)
        bideg=[0, 1], sing=[1], mult=1, factor=t+1/5*RootOf(_Z^2+5)
        bideg=[1, 1], sing=[1], mult=1, factor=t*s+3*t+s+1
        bideg=[1, 1], sing=[1], mult=1, factor=1/2*t+t*s-1/2*s
        bideg=[1, 0], sing=[1], mult=1, factor=s+2/3-1/3*RootOf(_Z^2+5)
        bideg=[1, 0], sing=[1], mult=1, factor=s+2/3+1/3*RootOf(_Z^2+5)
parametric type="D4b" value=[1/2, 0, 1, 0, 1/2, 1, 0, 0]
        bideg=[1, 0], sing=[1], mult=1, factor=s+1/13*RootOf(_Z^2+65)
        bideg=[1, 0], sing=[1], mult=1, factor=s-1/13*RootOf(_Z^2+65)
        bideg=[1, 1], sing=[1], mult=1, factor=t*s-5/6*t-1/3-1/6*s
        bideg=[0, 1], sing=[1], mult=1, factor=t+4/9+1/9*RootOf(_Z^2+65)
        bideg=[0, 1], sing=[1], mult=1, factor=t+4/9-1/9*RootOf(_Z^2+65)
        bideg=[1, 1], sing=[1], mult=1, factor=t*s+2*s+1
parametric type="D4c" value=[1/2, 1/2, 0, 0, 1/2, 0, 1/2, 0]
        bideg=[1, 0], sing=[1], mult=1, factor=2/3+RootOf(3*_Z^2+2*_Z+3)+s
        bideg=[1, 0], sing=[1], mult=1, factor=s-RootOf(3*_Z^2+2*_Z+3)
        bideg=[1, 1], sing=[1], mult=1, factor=t*s+3*t+s+1
        bideg=[0, 1], sing=[1], mult=1, factor=t-1/4-3/4*RootOf(3*_Z^2+2*_Z+3)
        bideg=[0, 1], sing=[1], mult=1, factor=t+1/4+3/4*RootOf(3*_Z^2+2*_Z+3)
        bideg=[1, 1], sing=[1], mult=1, factor=-t+t*s-2*s
parametric type="D4d" value=[1, -3/2, 0, 0, 1, 0, -3/2, 0]
        bideg=[0, 1], sing=[1], mult=1, factor=t+1/5*RootOf(_Z^2+145)
        bideg=[0, 1], sing=[1], mult=1, factor=t-1/5*RootOf(_Z^2+145)
        bideg=[1, 1], sing=[1], mult=1, factor=-2*t+t*s+2*s-5
        bideg=[1, 1], sing=[1], mult=1, factor=-1/3*t+t*s-7/3*s+1
        bideg=[1, 0], sing=[1], mult=1, factor=s-12/17-1/17*RootOf(_Z^2+145)
        bideg=[1, 0], sing=[1], mult=1, factor=s-12/17+1/17*RootOf(_Z^2+145)
parametric type="D5a" value=[1, 0, 0, 0, 1, 3/2, 0, 0]
        bideg=[0, 1], sing=[1], mult=1, factor=8/3+RootOf(3*_Z^2+8*_Z+3)+t
        bideg=[0, 1], sing=[1], mult=1, factor=t-RootOf(3*_Z^2+8*_Z+3)
        bideg=[0, 1], sing=[1], mult=1, factor=t-RootOf(17*_Z^2+24*_Z+17)
        bideg=[0, 1], sing=[1], mult=1, factor=24/17+RootOf(17*_Z^2+24*_Z+17)+t
        bideg=[1, 0], sing=[1], mult=2, factor=s-RootOf(1+_Z^2)
        bideg=[1, 0], sing=[1], mult=2, factor=s+RootOf(1+_Z^2)
parametric type="D5b" value=[1, 0, 0, 0, 1, 2, 0, 0]
        bideg=[0, 1], sing=[1], mult=1, factor=4/3+RootOf(3*_Z^2+4*_Z+3)+t
        bideg=[0, 1], sing=[1], mult=1, factor=t-RootOf(3*_Z^2+4*_Z+3)
        bideg=[1, 0], sing=[1], mult=2, factor=s-RootOf(1+_Z^2)
        bideg=[1, 0], sing=[1], mult=2, factor=s+RootOf(1+_Z^2)
parametric type="D5c" value=[1, 0, 0, 0, 1, 5/2, 0, 0]
        bideg=[0, 1], sing=[1], mult=1, factor=t-RootOf(5*_Z^2+8*_Z+5)
        bideg=[0, 1], sing=[1], mult=1, factor=8/5+RootOf(5*_Z^2+8*_Z+5)+t
        bideg=[0, 1], sing=[1], mult=1, factor=t+40/33+RootOf(33*_Z^2+40*_Z+33)
        bideg=[0, 1], sing=[1], mult=1, factor=t-RootOf(33*_Z^2+40*_Z+33)
        bideg=[1, 0], sing=[1], mult=2, factor=s+4/3+5/3*RootOf(5*_Z^2+8*_Z+5)
        bideg=[1, 0], sing=[1], mult=2, factor=s-4/3-5/3*RootOf(5*_Z^2+8*_Z+5)

memory used=1790.1MB, alloc=66.9MB, time=11.12
```


