#region License

// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn

#endregion

using System.Collections.Concurrent;
using System.Collections.Generic;

namespace ptmrs
{
    /// <summary>
    /// 
    /// </summary>
    public interface IDataConection
    {
        List<List<SpectraPackageItem>> GetNewDataPackage(int maxSizeOfPackage, out int numberOfSpectraPacked,
            out int numberOfPeptidesPacked);

        BlockingCollection<progressMessage> GetProgressMessageQueue();
    }
}