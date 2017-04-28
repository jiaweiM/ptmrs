#region License

// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn

#endregion

using System.Globalization;

namespace ptmrs
{
    public class Constants
    {
        public const int ID_NOT_SET = -1;

        public static readonly PTMSiteProbMode SITEPROBABILITY_MODE = PTMSiteProbMode.AtLeastNTime;

        public static readonly bool NO_MOD_ON_C_TERM_K_OR_R = false;

        public static CultureInfo EN_US_CULTURE_INFO = new CultureInfo("en-US");
    }
}