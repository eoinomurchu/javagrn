package grn.helper;

import grn.Grn;
import grn.Protein;
import grn.ProteinProducer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import java.util.ArrayList;

/**
 * A helper class providing output methods for GRNs
 */
public class GRNPrinter {

  /**
   *
   */
  public static void printState(Grn grn) {
    for (int i = 0; i < grn.tfProteins.length - grn.numberOfInputs; i++)
      System.out.print("TF"+i+" ");
    for (int i = 0; i < grn.numberOfInputs; i++)
      System.out.print("I"+i+" ");
    for (int i = 0; i < grn.pProteins.length; i++)
      System.out.print("P" + i + " ");
    System.out.print("\n");

    for (int i = 0; i < grn.tfProteins.length - grn.numberOfInputs; i++)
      System.out.print(grn.tfProteins[i].concentration+":"+i+" ");
    for (int i = grn.tfProteins.length - grn.numberOfInputs; i < grn.tfProteins.length; i++)
      System.out.print(grn.tfProteins[i].concentration+" ");
    for (int i = 0; i < grn.pProteins.length; i++)
      System.out.print(grn.pProteins[i].concentration+ " ");
    System.out.print("\n");
  }


  /**
   *
   */
  public static void printGRNToFile(final String fileName, final Grn grn, final double[][] results) {
    NumberFormat formatter = new DecimalFormat("###.#########");
    try{
      File file =new File(fileName);

      //if file doesnt exists, then create it
      boolean created = false;
      if(!file.exists()){
        file.createNewFile();
        created = true;
      }

      //true = append file
      FileWriter fileWritter = new FileWriter(file, true);
      BufferedWriter bufferWritter = new BufferedWriter(fileWritter);

      if (created) {
        for (int i = 0; i < grn.tfProteins.length - grn.numberOfInputs; i++)
          bufferWritter.write("TF"+i+" ");
        for (int i = 0; i < grn.numberOfInputs; i++)
          bufferWritter.write("I"+i+" ");
        for (int i = 0; i < grn.pProteins.length; i++)
          bufferWritter.write("P"+i+" ");
        bufferWritter.write("\n");
      }

      for (int t = 0; t < results.length; t++) {
        for (int i = 0; i < results[t].length; i++)
          bufferWritter.write(formatter.format(results[t][i])+" ");
        bufferWritter.write("\n");
      }
      bufferWritter.close();
      fileWritter.close();
    }
    catch(IOException e){
      e.printStackTrace();
    }
  }

  /**
   *
   */
  public static void printGRNasGraph(final Grn grn, final String fileName) {
    if (grn == null) {
      System.out.println("Trying to print null GRN");
      System.exit(-1);
    }

    int[][][] cbits = ProteinProducer.cbits;
    double[] exp    = ProteinProducer.exp;
    double edgeScale = 10.0;

    try{
      File file =new File(fileName+".dot");

      //if file doesnt exists, then create it
      if(!file.exists()){
        file.createNewFile();
      }

      //true = append file
      FileWriter fileWritter = new FileWriter(file, true);
      BufferedWriter bufferWritter = new BufferedWriter(fileWritter);

      bufferWritter.write("digraph graph_"+fileName+" {\n");

      bufferWritter.write("  // TF GENES\n");
      for (int i = 0; i < grn.tfGenes.length; i++) {
        for (int j = 0; j < grn.tfProteins.length - grn.numberOfInputs; j++) {
          double signal = exp[cbits[0][i][j]] - exp[cbits[1][i][j]];
          bufferWritter.write("  TF_"+j+" -> TF_"+i+" [color=\""+(signal < 0 ? "red" : "blue") + "\" penwidth=\""+Math.abs(signal*edgeScale)+"\"];\n");
        }
        for (int j = 0; j < grn.numberOfInputs; j++) {
          double signal = exp[cbits[0][i][grn.tfGenes.length+j]]-exp[cbits[1][i][grn.tfGenes.length+j]];
          bufferWritter.write("  I_"+j+" -> TF_"+i+" [color=\""+(signal < 0 ? "red" : "blue") + "\" penwidth=\""+Math.abs(signal*edgeScale)+"\"];\n");
        }
      }

      bufferWritter.write("  // P GENES\n");
      for (int i = 0; i < grn.pGenes.length; i++) {
        for (int j = 0; j < grn.tfProteins.length-grn.numberOfInputs; j++) {
          double signal = exp[cbits[0][grn.tfGenes.length+i][j]]-exp[cbits[1][grn.tfGenes.length+i][j]];
          bufferWritter.write("  TF_"+j+" -> P_"+i+" [color=\""+(signal < 0 ? "red" : "blue") + "\" penwidth=\""+Math.abs(signal*edgeScale)+"\"];\n");
        }
        for (int j = 0; j < grn.numberOfInputs; j++) {
          double signal = exp[cbits[0][grn.tfGenes.length+i][grn.tfGenes.length+j]]-exp[cbits[1][grn.tfGenes.length+i][grn.tfGenes.length+j]];
          bufferWritter.write("  I_"+j+" -> P_"+i+" [color=\""+(signal < 0 ? "red" : "blue") + "\" penwidth=\""+Math.abs(signal*edgeScale)+"\"];\n");
        }
      }

      bufferWritter.write("  { rank=source; ");
      for (int i = 0; i < grn.numberOfInputs; i++)
        bufferWritter.write("I_"+i+" ");
      bufferWritter.write("}\n");
      bufferWritter.write("  { rank=tfs; ");
      for (int i = 0; i < grn.tfProteins.length - grn.numberOfInputs; i++)
        bufferWritter.write("TF_"+i+" ");
      bufferWritter.write("}\n");
      bufferWritter.write("  { rank=sink; ");
      for (int i = 0; i < grn.pProteins.length; i++)
        bufferWritter.write("P_"+i+" ");
      bufferWritter.write("}\n");
      bufferWritter.write("}\n");
      bufferWritter.flush();
      bufferWritter.close();
      fileWritter.close();
    }
    catch(IOException e){
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    Grn grn = GRNLoader.readFromGenesFile(args[0]);
    ArrayList<Protein> ins = new ArrayList<Protein>();
    ins.add(new Protein(0.0, 0x00000000));
    ins.add(new Protein(0.1, 0x0000FFFF));
    ins.add(new Protein(0.1, 0xFFFF0000));
    ins.add(new Protein(0.1, 0xFFFFFFFF));
    grn.injectInputs(ins);
    GRNPrinter.printGRNasGraph(grn, "graph");
  }
}
