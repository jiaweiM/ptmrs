#region License

// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn

#endregion

using System;

namespace ptmrs
{
    public class StringValue : Attribute
    {
        public string Value { get; set; }

        public StringValue(string str)
        {
            Value = str;
        }
    }
}