#region License
// Created at 2016.07.13
// Author:  JiaweiMao
// Copyright (c) Dalian Institute of Chemical Physics
//            Chinese Academy of Sciences
// Contact: jiawei@dicp.ac.cn
#endregion
namespace ptmrs
{
    /// <summary>
    /// 
    /// </summary>
    public class progressMessage
    {
        public enum typeOfMessage
        {
            noType,
            progressMessage,
            stringMessage
        }

        public double spectraProcessed;

        public int numberOfPeptidesProcessed;

        public string message;

        private typeOfMessage msgType;

        public typeOfMessage type
        {
            get
            {
                return this.msgType;
            }
            private set
            {
                this.msgType = value;
            }
        }

        public progressMessage(string message)
        {
            this.message = message;
            this.spectraProcessed = 0.0;
            this.numberOfPeptidesProcessed = 0;
            this.msgType = progressMessage.typeOfMessage.stringMessage;
        }

        public progressMessage(double spectra, int peptide)
        {
            this.message = null;
            this.spectraProcessed = spectra;
            this.numberOfPeptidesProcessed = peptide;
            this.msgType = progressMessage.typeOfMessage.progressMessage;
        }
    }
}