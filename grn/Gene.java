package grn;

/**
 *
 */
public class Gene {

  /** The promoter signature for TF genes */
  public static int TF_PROMOTER = 0x00000000;

  /** The promoter signature for P genes */
  public static int P_PROMOTER = 0x000000FF;

  /** The size of a gene in 32 bit units */
  public static int SIZE = 8;

  /** The gene's enhancer regulatory site signature */
  public int enhancer;

  /** The gene's inhibitor regulatory site signature */
  public int inhibitor;

  /** The gene's promoter site signature */
  public int promoter;

  /** The codons making up the gene's protein coding region */
  public int codons[];

  /** An index field */
  public int index;

  /**
   * Default constructor - not implemented
   */
  public Gene() {
    ;
  }

}
