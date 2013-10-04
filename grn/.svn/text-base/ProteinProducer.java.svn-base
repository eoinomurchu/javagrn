package grn;

import java.util.ArrayList;
import java.util.HashMap;

public class ProteinProducer {

  public static double[] exp;

  /** */
  private static double scalingFactor = 1.0;

  /** */
  private static double delta = 1.0;

  /** */
  private static double phi(double d) {
    return 0.0;
  }

  private static HashMap<Integer, Integer> compBitsMap = new HashMap<Integer, Integer>();

  public static int[][][] cbits;

  static {
    exp = new double[33];
    for (int i = 0; i < 33; i++)
      exp[i] = Math.exp(-i);
  }

  /**
   *
   *
   * @param genes
   * @return Proteins
   */
  public static Protein[] expressGenes(Gene[] genes) {
    Protein[] proteins = new Protein[genes.length];
    int nP = 0;

    for (Gene gene : genes) {
      Protein p = new Protein();

      int[] bits = new int[32];
      for (int i = 0; i < 5; i++) {
        int codon = gene.codons[i];
        for (int j = 0; j < 32; j++) {
          bits[j] += (codon >>> j) & 1;
        }
      }

      p.value = 0;
      for (int j = 31; j > -1; j--) {
        p.value <<= 1;
        p.value += bits[j] > 2 ? 1 : 0;
      }
      proteins[nP++] = p;
    }
    return proteins;
  }

  /**
   *
   *
   * @param g
   * @param p
   * @param proteins
   * @return production delta
   */
  public static double produce(Gene g, Protein p, Protein[] proteins, int umax) {
    return (delta * (enhancerSignal(g, proteins, umax) - inhibitorSignal(g, proteins, umax)) * p.concentration - phi(1.0));
  }

  public static double p_produce(Gene g, Protein p, Protein[] proteins, int umax) {
    return (delta * (enhancerSignal(g, proteins, umax) - inhibitorSignal(g, proteins, umax)) - phi(1.0));
  }


  /**
   *
   *
   * @param
   * @param
   * @return
   */
  private static double enhancerSignal(Gene g, Protein[] proteins, int umax) {
    return regulatorySignal(g.enhancer, proteins, umax, g.index, 0);
  }

  /**
   *
   *
   * @param
   * @param
   * @return
   */
  private static double inhibitorSignal(Gene g, Protein[] proteins, int umax) {
    return regulatorySignal(g.inhibitor, proteins, umax, g.index, 1);
  }

  /**
   *
   *
   * @param
   * @param
   * @return
   */
  private static double regulatorySignal(int r, Protein[] proteins, int umax, int index, int type) {
    double signal = 0.0;

    // int[] cbits = new int[proteins.length];
    // for (int i = 0; i < proteins.length; i++) {
    //   cbits[i] = countComplementaryBits(r, proteins[i]);
    // }

    for (int i = 0; i < proteins.length; i++) {

       int compMeasure;
       if (cbits[type][index][i] == -1) {
         compMeasure = umax - countComplementaryBits(r, proteins[i]);
         cbits[type][index][i] = compMeasure;
       }
       else
         compMeasure = cbits[type][index][i];

       signal += proteins[i].concentration * exp[compMeasure];
       //signal += proteins[i].concentration * exp[umax-cbits[i]];//Math.exp(scalingFactor * (cbits[i] - max_cbits));
    }

    if (proteins.length == 0)
      return 0.0;

    return signal/(double)proteins.length;
  }

  /**
   *
   *
   * @param
   * @param
   * @return
   */
  private static int countComplementaryBits(int r, Protein p) {

    int c = 0;
    int v = r ^ p.value;
    // while (v != 0) {
    //   c++;
    //   v &= v - 1;
    // }
    v = v - ((v >>> 1) & 0x55555555);                    // reuse input as temporary
    v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);     // temp
    c = ((v + (v >>> 4) & 0xF0F0F0F) * 0x1010101) >>> 24; // count
    return c;
  }

  public static int calculateUMax(Gene g, Protein[] proteins) {
    int max_cbits = -1, temp;
    for (int i = 0; i < proteins.length; i++) {
      temp = countComplementaryBits(g.enhancer, proteins[i]);
      if(temp > max_cbits)
        max_cbits = temp;
      temp = countComplementaryBits(g.inhibitor, proteins[i]);
      if(temp > max_cbits)
        max_cbits = temp;
    }
    return max_cbits;
  }
}
