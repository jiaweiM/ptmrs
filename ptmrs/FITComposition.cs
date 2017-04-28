using System;

namespace ptmrs
{
    [Flags]
    public enum FitComposition
    {
        No = 0,
        B = 1,
        Y = 2,
        C = 4,
        ZRadical = 8,
        ZPrime = 16
    }
}