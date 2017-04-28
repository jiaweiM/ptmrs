#region License

// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn

#endregion

using System.Collections.Generic;

namespace ptmrs
{
    /// <summary>
    /// 
    /// </summary>
    internal static class FragmentTypeExtensionMethods
    {
        public static int Contains(this FragmentIonType type, bool ignorePriority, bool ignoreExpansion,
            bool ignoreAppendix, List<FragmentIonType> set)
        {
            for (int i = 0; i < set.Count; i++)
            {
                if (set[i] != null && set[i].Equals(type, ignorePriority, ignoreExpansion, ignoreAppendix, false, false))
                {
                    return i;
                }
            }
            return -1;
        }

        public static int TryGetFITIndex(this List<FragmentIonType> set, FragmentIonType value)
        {
            for (int i = 0; i < set.Count; i++)
            {
                if (set[i].Equals(value, true, true, true, false, true))
                {
                    return i;
                }
            }
            return -1;
        }
    }
}