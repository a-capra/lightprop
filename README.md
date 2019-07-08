# lightprop

I need [ROOT][rootlink]

[rootlink]: https://root.cern.ch/



Compile me like this:

```
g++ -Wall -O3 `root-config --cflags` `root-config --glibs` -o LightProp LightProp.cc
```