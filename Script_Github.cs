using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.IO;
using System.Data;

using ILOG.Concert;
using ILOG.CPLEX;

namespace Example1
{
    public class Data
    {
        public int[][] mutMatrix = null;
        public double[][] connectMat = null;

        // Constructor.
        public Data()
        {
        }

        public void ReadData(string arquivo)
        {
            int lin, col = 0;

            // Opens file to read data.
            StreamReader ObjReader = new StreamReader(arquivo, Encoding.UTF7);

            string sLine = "";

            // Reads dimension of mutation matrix.
            sLine = ObjReader.ReadLine();

            string[] variavel = sLine.Split('\t');

            lin = Convert.ToInt32(variavel[0]);
            col = Convert.ToInt32(variavel[1]);

            // Allocates memory for mutation matrix.
            this.mutMatrix = new int[lin][];
            for (int i = 0; i < lin; i++)
                this.mutMatrix[i] = new int[col];

            // Loads entries of mutation matrix.
            for (int i = 0; i < lin; i++)
            {
                sLine = ObjReader.ReadLine();
                variavel = sLine.Split('\t');

                for (int j = 0; j < col; j++)
                    this.mutMatrix[i][j] = Convert.ToInt32(variavel[j]);
            }

            // Skips one line.
            sLine = ObjReader.ReadLine();

            // Reads dimension of connectivity matrix.
            sLine = ObjReader.ReadLine();

            variavel = sLine.Split('\t');

            lin = Convert.ToInt32(variavel[0]);
            col = Convert.ToInt32(variavel[1]);

            // Allocates memory for connectivity matrix.
            this.connectMat = new double[lin][];
            for (int i = 0; i < lin; i++)
                this.connectMat[i] = new double[col];

            // Loads entries of connectivity matrix.
            for (int i = 0; i < lin; i++)
            {
                sLine = ObjReader.ReadLine();
                variavel = sLine.Split('\t');

                for (int j = 0; j < col; j++)
                    this.connectMat[i][j] = Convert.ToDouble(variavel[j]);
            }

            // Closes file.
            ObjReader.Close();
        }
    }

    class Program
    {
        static void Main(string[] args)
        {

            int nSamples = 0; // Number of samples
            int nMutGenes = 0; // Number of mutation genes
            int nExpGenes = 0; // Number of expression genes
            int nVar = 0; // Number of variables in MILP
            int nPathways = 5; // Number of pathways
            int nConstraints = 0; // Number of constraints in MILP

            Data data = new Data();

            // Loads input data.
            data.ReadData("../../inputFile.txt");

            nSamples = data.mutMatrix.Length;
            nMutGenes = data.mutMatrix[0].Length;
            nExpGenes = data.connectMat.Length;

            // Total number of decision variables:
            // pM --> nMutGenes * nPathways
            // pE --> nExpGenes * nPathways
            // aM --> nSamples * nPathways
            // fM --> nSamples * nPathways
            nVar = (nMutGenes * nPathways) + (nExpGenes * nPathways) + 2 * (nSamples * nPathways);

            // Instance of cplex model.
            Cplex cplex = new Cplex();

            /* We were having issues with "Out of Memory" errors as we increased the number of samples and number of genes
             * To address this, we set parameter NodeFileInd = 3 --> this indicates that node info from MIP optimizer will be compressed 
             * and written to files instead of accessing the computer memory during the simulation run
             *
             * */
            cplex.SetParam(Cplex.IntParam.NodeFileInd, 3);

            // Sets names of decision variables.
            string[] name = new string[nVar];

            for (int i = 1; i <= nVar; i++)
            {
                name[i - 1] = "x" + i.ToString();
            }

            // Sets type, lower and upper bound of decision variables.
            NumVarType[] varType = new NumVarType[nVar];
            double[] lb = new double[nVar];
            double[] ub = new double[nVar];

            for (int i = 0; i < nMutGenes * nPathways; i++)
            {
                // Variable pM is binary.
                varType[i] = NumVarType.Bool;
                lb[i] = 0;
                ub[i] = 1;
            }

            for (int i = nMutGenes * nPathways; i < (nMutGenes * nPathways) + (nExpGenes * nPathways); i++)
            {
                // Variable pE is real number.
                varType[i] = NumVarType.Float;
                lb[i] = 0;
                ub[i] = System.Double.MaxValue;
            }

            for (int i = (nMutGenes * nPathways) + (nExpGenes * nPathways); i < nVar; i++)
            {
                // Variables aM and fM are binary.
                varType[i] = NumVarType.Bool;
                lb[i] = 0;
                ub[i] = 1;
            }


            // Decision Variables.
            INumVar[] x = cplex.NumVarArray(nVar, lb, ub, varType, name);

            // Coeficients of objective function.
            double[] objvals = new double[nVar];

            int contobjvals = 0;
            int idMap = 0; // Mapping index.

            // Variable pM.
            int[][] pM = new int[nMutGenes][]; // Mapping matrix pM.
            for (int i = 0; i < nMutGenes; i++)
            {
                pM[i] = new int[nPathways];
                for (int j = 0; j < nPathways; j++)
                {
                    // Mapping....
                    pM[i][j] = idMap;
                    idMap++;

                    objvals[contobjvals] = 0;
                    for (int k = 0; k < nSamples; k++)
                        // Weight associated with first component of objective function = 0.1.
                        objvals[contobjvals] += 0.1 * data.mutMatrix[k][i];

                    contobjvals++;
                }
            }

            // Variable pE.
            int[][] pE = new int[nExpGenes][];
            for (int i = 0; i < nExpGenes; i++)
            {
                pE[i] = new int[nPathways];
                for (int j = 0; j < nPathways; j++)
                {
                    // Mapping....
                    pE[i][j] = idMap;
                    idMap++;

                    // Weight associated with second component of objective function = 0.9.
                    objvals[contobjvals] = -0.9;
                    contobjvals++;
                }
            }

            // Variable aM.
            int[][] aM = new int[nSamples][];
            for (int i = 0; i < nSamples; i++)
            {
                aM[i] = new int[nPathways];
                for (int j = 0; j < nPathways; j++)
                {
                    // Mapping....
                    aM[i][j] = idMap;
                    idMap++;

                    // Weight associated with first component of objective function = 0.1.
                    objvals[contobjvals] = -0.1;
                    contobjvals++;
                }
            }

            // Variable fM.
            int[][] fM = new int[nSamples][];
            for (int i = 0; i < nSamples; i++)
            {
                fM[i] = new int[nPathways];
                for (int j = 0; j < nPathways; j++)
                {
                    // Mapping....
                    fM[i][j] = idMap;
                    idMap++;

                    // Weight associated with first component of objective function = 0.1.
                    objvals[contobjvals] = 0.1;
                    contobjvals++;
                }
            }

            // Defines the cost function as a minimization of decision variables with corresponding cost coefficients.
            cplex.Add(cplex.Minimize(cplex.ScalProd(x, objvals)));

            nConstraints = nMutGenes + nExpGenes + 2 * nPathways + nSamples * (nPathways - 1) + (nSamples * nPathways) + (nPathways * nExpGenes); // Total number of constraints.
            IRange[] constraint = new IRange[nConstraints];
            int contconstraints = 0;
            INumExpr[] Expr = null;

            // Constraints C1.
            Expr = new INumExpr[nPathways];
            for (int i = 0; i < nMutGenes; i++)
            {
                for (int j = 0; j < nPathways; j++)
                    Expr[j] = x[pM[i][j]];

                // Constraints of the form sum(pM) = 1.
                constraint[contconstraints] = cplex.AddRange(1, 1); // Defines right hand side (= 1).
                constraint[contconstraints].Expr = cplex.Sum(Expr); // Defines left hand side (sum(pM)).
                contconstraints++;
            }

            // Constraints C2.
            Expr = new INumExpr[nPathways];
            for (int i = 0; i < nExpGenes; i++)
            {
                for (int j = 0; j < nPathways; j++)
                    Expr[j] = x[pE[i][j]];

                // Constraints of the form sum(pE) > 0.
                constraint[contconstraints] = cplex.AddRange(0.000001, System.Double.MaxValue); // Defines right hand side (> 0).
                constraint[contconstraints].Expr = cplex.Sum(Expr); // Defines left hand side (sum(pE)).
                contconstraints++;
            }

            // Constraints C3.
            Expr = new INumExpr[nMutGenes];
            for (int i = 0; i < nPathways; i++)
            {
                for (int j = 0; j < nMutGenes; j++)
                    Expr[j] = x[pM[j][i]];

                // Constraints of the form sum(pM) >= 1.
                constraint[contconstraints] = cplex.AddRange(1, System.Double.MaxValue); // Defines right hand side (>= 1).
                constraint[contconstraints].Expr = cplex.Sum(Expr); // Defines left hand side (sum(pM)).
                contconstraints++;
            }

            // Constraints C4.
            int id = nMutGenes * nPathways;
            Expr = new INumExpr[nExpGenes];
            for (int i = 0; i < nPathways; i++)
            {
                for (int j = 0; j < nExpGenes; j++)
                    Expr[j] = x[pE[j][i]];

                // Constraints of the form sum(pE) > 0.
                constraint[contconstraints] = cplex.AddRange(0.000001, System.Double.MaxValue); // Defines right hand side (> 0).
                constraint[contconstraints].Expr = cplex.Sum(Expr); // Defines left hand side (sum(pE)).
                contconstraints++;
            }

            // Constraints C5.
            id = nMutGenes * nPathways + nExpGenes * nPathways;
            for (int i = 0; i < nSamples; i++)
            {
                // Constraints of the form aM - aM >= 0.
                for (int j = 0; j < (nPathways - 1); j++)
                {
                    constraint[contconstraints] = cplex.AddRange(0, System.Double.MaxValue); // Defines right hand side (>= 0).
                    constraint[contconstraints].Expr = cplex.Sum(cplex.Prod(1.0, x[aM[i][j]]), cplex.Prod(-1.0, x[aM[i][j + 1]])); // Defines left hand side (aM - aM).
                    contconstraints++;
                }
            }

            // Constraints C6.
            Expr = new INumExpr[nMutGenes + 2];
            for (int i = 0; i < nSamples; i++) // Goes along lines of mutation matrix.
            {
                for (int k = 0; k < nPathways; k++)
                {
                    for (int j = 0; j < nMutGenes; j++) // Goes along columns of mutation matrix.
                    {
                        Expr[j] = cplex.Prod(data.mutMatrix[i][j], x[pM[j][k]]);
                    }

                    Expr[nMutGenes] = cplex.Prod(1.0, x[fM[i][k]]);
                    Expr[nMutGenes + 1] = cplex.Prod(-1.0, x[aM[i][k]]);

                    // Constraints of the form sum(mutMat*pM) + fM - aM >= 0.
                    constraint[contconstraints] = cplex.AddRange(0, System.Double.MaxValue); // Defines right hand side (>= 0).
                    constraint[contconstraints].Expr = cplex.Sum(Expr); // Defines left hand side (sum(mutMat*pM) + fM - aM).
                    contconstraints++;
                }
            }

            // Constraints C7.
            Expr = new INumExpr[nMutGenes + 1];
            for (int i = 0; i < nPathways; i++)
            {
                for (int j = 0; j < nExpGenes; j++)
                {
                    for (int k = 0; k < nMutGenes; k++)
                    {
                        Expr[k] = cplex.Prod(-data.connectMat[j][k], x[pM[k][i]]);
                    }

                    Expr[nMutGenes] = cplex.Prod(1.0, x[pE[j][i]]);
                    
                    // Constraints of the form pE - sum(connectMat*pM) = 0.
                    constraint[contconstraints] = cplex.AddRange(0, 0); // Defines right hand side (= 0).
                    constraint[contconstraints].Expr = cplex.Sum(Expr); // Defines left hand side (pE - sum(connectMat*pM)).
                    contconstraints++;
                }
            }

            // Writes file with MILP model formulation.
            cplex.ExportModel("MILP.lp");

            try
            {
                if (cplex.Solve())
                {
                    // Retrieves optimal values of decision variables.
                    double[] xx = cplex.GetValues(x);

                    // Displays solution info on command window.
                    cplex.Output().WriteLine("Solution status=" + cplex.GetStatus());
                    cplex.Output().WriteLine("Solution value = " + cplex.ObjValue);

                    // Writes output to txt file.
                    StreamWriter objWriter = new StreamWriter("Result.txt");
                    int cont = 0;
                    string sLine = "";
                    int nvars = xx.Length;

                    objWriter.WriteLine("Total cost");
                    objWriter.WriteLine(cplex.ObjValue.ToString());
                    objWriter.WriteLine("");

                    objWriter.WriteLine("Matrix pM");
                    objWriter.WriteLine("");
                    for (int i = 0; i < nMutGenes; i++)
                    {
                        sLine = "";
                        for (int j = 0; j < nPathways; j++)
                        {
                            sLine += xx[cont].ToString() + "\t";
                            cont++;
                        }

                        objWriter.WriteLine(sLine);
                    }

                    // Value of second term of objective function (sum of pE's).
                    double sumPE = 0.0;

                    objWriter.WriteLine("");
                    objWriter.WriteLine("Matrix pE");
                    objWriter.WriteLine("");
                    for (int i = 0; i < nExpGenes; i++)
                    {
                        sLine = "";
                        for (int j = 0; j < nPathways; j++)
                        {
                            sLine += xx[cont].ToString() + "\t";
                            cont++;

                            // Computes running sum of pE's.
                            sumPE = sumPE + xx[cont];
                        }

                        objWriter.WriteLine(sLine);
                    }

                    objWriter.WriteLine("");
                    objWriter.WriteLine("Sum pE");
                    objWriter.WriteLine(sumPE.ToString());
                    objWriter.WriteLine("");

                    objWriter.WriteLine("");
                    objWriter.WriteLine("Matrix aM");
                    objWriter.WriteLine("");
                    for (int i = 0; i < nSamples; i++)
                    {
                        sLine = "";
                        for (int j = 0; j < nPathways; j++)
                        {
                            sLine += xx[cont].ToString() + "\t";
                            cont++;
                        }

                        objWriter.WriteLine(sLine);
                    }

                    objWriter.WriteLine("");
                    objWriter.WriteLine("Matrix fM");
                    objWriter.WriteLine("");
                    for (int i = 0; i < nSamples; i++)
                    {
                        sLine = "";
                        for (int j = 0; j < nPathways; j++)
                        {
                            sLine += xx[cont].ToString() + "\t";
                            cont++;
                        }

                        objWriter.WriteLine(sLine);
                    }

                    objWriter.Close();

                    // Writes optimal solution info to file.
                    cplex.WriteSolution("solution");
                }
                else 
                {
                    cplex.GetCplexStatus();
                }
                cplex.End();
            }
            catch (ILOG.Concert.Exception e)
            {
                System.Console.WriteLine("Concert exception '" + e + "' caught");
            }
        }
    }
}