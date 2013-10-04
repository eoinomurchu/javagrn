package grn;

import java.util.ArrayList;

import grn.helper.ArrayUtils;

/**
 *
 */
public class BitScanner {

  /** A bit mask for a promoter site (the last 8 bits) */
  public static int PROMO_MASK = 0x000000FF;

  /** The bit index to start searching for genes from  */
  private static int STARTING_INDEX = -(Gene.SIZE - 2)*32;

  /** The current search index */
  private int currentIndex;

  /** The chromosome/genome being searched */
  public int[] chromo;

  /** The length of the chromosome being searched */
  public int chromoLength;

  /** The TF genes located */
  Gene[] tfGenes;

  /** The number of TF genes located */
  int nTfGenes;

  /** The P genes located */
  Gene[] pGenes;

  /** The number of P genes located */
  int nPGenes;

  /**
   * Creates a new BitScanner and initiates the search for genes.
   * 
   * @param codons the genome to be searched - a bit sequence represented by 32 bit integers
   */
  public BitScanner(int[] codons) {
    chromo = codons;
    chromoLength = chromo.length;

    tfGenes = new Gene[chromoLength/8];
    pGenes = new Gene[chromoLength/8];

    currentIndex = STARTING_INDEX;

    findGenes();
  }

  /**
   * @return an array of TF genes found along the genome
   */
  public Gene[] getTFGenes() {
    return tfGenes;
  }

  /**
   * @return an array of P genes found along the genome
   */
  public Gene[] getPGenes() {
    return pGenes;
  }

  /**
   * Performs the search for genes along the chromosome
   */
  private void findGenes() {
    nTfGenes = 0;
    nPGenes  = 0;

    /* Process genes while nextPromoter() returns a positive value */
    while ((currentIndex = nextPromoter(currentIndex)) > -1) {
      /* Check if there is sufficient storage for more genes */
      tfGenes = ArrayUtils.checkAndResizeArray(tfGenes, nTfGenes);
      pGenes = ArrayUtils.checkAndResizeArray(pGenes, nTfGenes);

      /* Check if the current bit is the beginning of a promoter */
      if (isAPromoter(currentIndex, Gene.TF_PROMOTER)){
        tfGenes[nTfGenes++] = getGeneFromPromoterIndex(currentIndex);
      }
      else if (isAPromoter(currentIndex, Gene.P_PROMOTER)){
        pGenes[nPGenes++] = getGeneFromPromoterIndex(currentIndex);
      }
    }

    /* Resize over sized arrays */
    tfGenes = ArrayUtils.fitArray(tfGenes, nTfGenes);
    pGenes  = ArrayUtils.fitArray(pGenes, nPGenes);
  }

  /**
   * Constructs a 32 bit int from 32 bits along the bit sequence from
   * the current index. The bit sequence is made up using 32 bit
   * integers. The returned integer can contain bits from two
   * neighbouring integers.
   *
   * @param index a bit index along the bit string genome
   * @return the int value representing the 32 bits from index
   */
  private int getIntFromBitIndex(int index) {
    /* Calculate which codon and how far along it the index is */
    int codon = index / 32;
    int bitIndex = index % 32;

    /* 
     * Index is at the start of an integer so just return that int.
     * Special case required as >>> is mod 32
     */
    if (bitIndex == 0)
      return chromo[codon];
    /* Otherwise just construct the int from segments of this codon and the next codon */
    else if (codon < chromoLength - 1) {
      return (chromo[codon] << bitIndex) | (chromo[codon+1] >>> (32 - bitIndex));
    }
    else {
      System.out.println("Error getting codon value, index too near end of data: "+index+" -> "+codon+":"+bitIndex+" > "+(chromoLength -1));
      return -1;
    }
  }

  /**
   * Get the bit at a specific bit index along the bit string genome
   * 
   * @param index bit index
   * @return the value of the bit, 0 or 1
   */
  private int getBit(int index) {
    int codon = index / 32;
    int bitIndex = index % 32;

    int bit = chromo[codon] << bitIndex;
    return bit >>> 31;
  }

  /**
   * Construct a gene from a promoter's bit index
   *
   * @param index beginning of promoter 
   * @return The gene about the promoter index
   */
  private Gene getGeneFromPromoterIndex(int index) {
    Gene g = new Gene();
    g.enhancer  = getIntFromBitIndex(index-64);
    g.inhibitor = getIntFromBitIndex(index-32);
    g.promoter  = getIntFromBitIndex(index);
    g.codons    = new int[5];

    for (int i = 0; i < 5; i++)
      g.codons[i] = getIntFromBitIndex(index + 32 + (i * 32));

    return g;
  }

  /**
   * Find the next promoter site of either gene class. Overlapping
   * genes are omitted.
   *
   * @param prevIndex the index of the previous promoter site 
   * @return the index of the next promoter site
   */
  private int nextPromoter(int prevIndex) {
    if (prevIndex < 0)
      prevIndex = STARTING_INDEX;

    /* Jump over the previous gene */
    prevIndex += Gene.SIZE*32;

    /* Increment along the bit sequence while no promoter has been
     * found, and there are sufficient bits remaining */
    while (prevIndex <= chromoLength*32 - 192 &&
           !isAPromoter(prevIndex, Gene.TF_PROMOTER) &&
           !isAPromoter(prevIndex, Gene.P_PROMOTER))
      prevIndex++;

    /* Not sufficient bits left to create anymore genes */
    if (prevIndex > chromoLength*32 - 192)
      prevIndex = -1;

    return prevIndex;
  }


  /**
   * Check if there is a promoter at the current index. The 32 bits
   * starting from index are gathered and used to check against the
   * promoter mask.
   *
   * @param index the bit index
   * @param mask a gene class promoter mask
   * @return true if a promoter is detected
   */
  private boolean isAPromoter(int index, int mask) {
    int i = getIntFromBitIndex(index);
    return ((i & PROMO_MASK) ^ mask) == 0;
  }
}
