#region License

// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn

#endregion

using System.ComponentModel;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Resources;

namespace ptmrs.Properties
{
    public class Resources
    {
        private static ResourceManager resourceMan;
        private static CultureInfo resourceCulture;

        [EditorBrowsable(EditorBrowsableState.Advanced)]
        public static ResourceManager ResourceManager
        {
            get
            {
                if (object.ReferenceEquals(Resources.resourceMan, null))
                {
                    ResourceManager resourceManager = new ResourceManager("IMP.ptmRS.Properties.Resources",
                        typeof (Resources).Assembly);
                    Resources.resourceMan = resourceManager;
                }
                return Resources.resourceMan;
            }
        }

        [EditorBrowsable(EditorBrowsableState.Advanced)]
        public static CultureInfo Culture
        {
            get { return Resources.resourceCulture; }
            set { Resources.resourceCulture = value; }
        }

        public static string IMP_ptmRSConf
        {
            get { return Resources.ResourceManager.GetString("IMP_ptmRSConf", Resources.resourceCulture); }
        }

        [SuppressMessage("Microsoft.Performance", "CA1811:AvoidUncalledPrivateCode")]
        internal Resources() {}
    }
}