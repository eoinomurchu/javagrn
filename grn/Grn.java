package grn;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import grn.helper.ArrayUtils;
import grn.helper.GRNPrinter;

/**
 *
 */
public class Grn {

  /** This is the lowest limit of concentration for a protein */
  public static final double ZERO = 1e-10;

  /** The length of the initial stablisation syncronisation step */
  private static final int INIT_PERIOD = 10000;

  /** Widow size when checking for stability/steady state */
  private static final int REST_STEP = 100;

  /** Threshold of stability */
  private static final double REST_EPSILON = ZERO;

  /** 
   * The GRN time series data for the latest syncronisation step.
   * Protein concentration values are stored at each time step:
   * t x (TF proteins, Input TF proteins, P proteins)
   */
  public double[][] results;

  /** The GRN time series data for stabilisation */
  public double[][] initResults;

  /** The genome of integers */
  public int[] chromosome;

  /** Transcription Factor Genes */
  public Gene[] tfGenes;

  /** Product Genes */
  public Gene[] pGenes;

  /** Transcription Factor Proteins */
  public Protein[] tfProteins;

  /** Product Proteins */
  public Protein[] pProteins;

  /** Sum total of input TF Protein concentrations */
  public double inputConcentration;

  /** Number of input TF Proteins */
  public int numberOfInputs;

  /** The greatest match observed between a protein and a regulatory site */
  public int umax;

  /**
   * Creates a new GRN from a binary string in the form of an integer
   * array. No inputs proteins are provided.
   *
   * @param codons 32 bit integers making up the binary string
   */
  public Grn(int[] codons) {
    this(codons, new ArrayList<Protein>());
  }

  /**
   * Creates a new GRN from a binary string in the form of an integer
   * array. Input proteins are provided and taken into account when
   * setting initial concentration values.
   *
   * The binary string is scanned, locating and extracting TF genes
   * and P genes. These genes are then expressed, producing sets of TF
   * proteins and P proteins.
   *
   * Initial concentration values are assigned to the proteins, and
   * the input proteins are added to the list of TF proteins.
   *
   * Finally the degree of matching between each TF protein and
   * regulatory site is precalculated.
   *
   * @param codons 32 bit integers making up the binary string
   * @param inputProteins list of input TF proteins
   */
  public Grn(int[] codons, ArrayList<Protein> inputProteins) {
    /* Locate Genes */
    BitScanner hunter = new BitScanner(codons);
    tfGenes = hunter.getTFGenes();
    pGenes = hunter.getPGenes();
    hunter = null;

    /* Express Proteins */
    tfProteins = ProteinProducer.expressGenes(tfGenes);
    pProteins = ProteinProducer.expressGenes(pGenes);

    /* Initialise concentrations */
    calculateInputConcentration(inputProteins);
    setInitialProteinConcentrations();

    /* Add inputs */
    numberOfInputs = inputProteins.size();
    tfProteins = Arrays.copyOf(tfProteins, tfProteins.length + numberOfInputs);
    for (int p = 0; p < numberOfInputs; p++)
      tfProteins[tfGenes.length + p] = inputProteins.get(p);

    /* Generate the precalculated tables */
    generateTables();
  }

  /**
   * Creates a new GRN from a binary string in the form of an integer
   * string. Input proteins are provided and taken into account when
   * setting initial concentration values.
   *
   * The binary string is scanned, locating and extracting TF genes
   * and P genes. These genes are then expressed, producing sets of TF
   * proteins and P proteins.
   *
   * Initial concentration values are assigned to the proteins, and
   * the input proteins are added to the list of TF proteins.
   *
   * Finally the degree of matching between each TF protein and
   * regulatory site is precalculated.
   *
   * @param codonString string of integer values encoding the binary genome
   * @param inputProteins list of input TF proteins
   */
  public Grn(String codonString, ArrayList<Protein> inputProteins) {
    /* Decode the string into ints for locating the genes */
    String[] codonStrings = codonString.split(" ");
    int[] codons = new int[codonStrings.length];
    for (int i = 0; i < codonStrings.length; i++)
      codons[i] = Integer.parseInt(codonStrings[i]);

    /* Locate Genes */
    BitScanner hunter = new BitScanner(codons);
    tfGenes = hunter.getTFGenes();
    pGenes = hunter.getPGenes();
    hunter = null;

    /* Express Proteins */
    tfProteins = ProteinProducer.expressGenes(tfGenes);
    pProteins = ProteinProducer.expressGenes(pGenes);

    /* Initialise concentrations */
    calculateInputConcentration(inputProteins);
    setInitialProteinConcentrations();

    /* Add inputs */
    numberOfInputs = inputProteins.size();
    tfProteins = Arrays.copyOf(tfProteins, tfProteins.length + numberOfInputs);
    for (int p = 0; p < numberOfInputs; p++)
      tfProteins[tfGenes.length + p] = inputProteins.get(p);

    /* Generate the precalculated tables */
    generateTables();
  }

  /**
   * Precalculate a lookup table of the degree of matching between
   * each TF protein and regulatory site in the network. 
   *
   * 3D lookup table: TF proteins x all genes x two regulatory sites per gene.
   *
   * Matching is calculated as the number of complementary bits
   * between a protein's signature and a regulatory site's value. This
   * is in the range of [0,32].
   *
   * The lookup table is populated with umax - complementaryBitCount.
   */
  private void generateTables() {
    /* Initialise the cbits array */
    ProteinProducer.cbits = new int[2][tfGenes.length+pGenes.length][tfProteins.length];
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < tfGenes.length + pGenes.length; j++)
        for (int k = 0; k < tfProteins.length; k++)
          ProteinProducer.cbits[i][j][k] = -1;

    /* Find the maximum level of matching */
    umax = -1;
    for (Gene g : tfGenes) {
      int temp;
      if ((temp = ProteinProducer.calculateUMax(g, tfProteins)) > umax)
        umax = temp;
    }
    for (Gene g : pGenes) {
      int temp;
      if ((temp = ProteinProducer.calculateUMax(g, tfProteins)) > umax)
        umax = temp;
    }

    /* Populate the table by expressing each gene once */
    for (int i = 0; i < tfGenes.length;  i++) {
      tfGenes[i].index = i;
      ProteinProducer.produce(tfGenes[i], tfProteins[i], tfProteins, umax);
    }
    for (int i = 0; i < pGenes.length;  i++) {
      pGenes[i].index = tfGenes.length+i;
      ProteinProducer.p_produce(pGenes[i], pProteins[i], tfProteins, umax);
    }
  }

  /**
   * Calculate the total input concentration
   */
  private void calculateInputConcentration(ArrayList<Protein> inputProteins) {
    inputConcentration = 0.0;
    for (Protein p : inputProteins)
      inputConcentration += p.concentration;
  }

  /**
   * Assign initial protein concentration values. Value is 1/N, where
   * N is the number of proteins in a class.  If there are inputs then
   * (1-inputs)/N is used.
   */
  private void setInitialProteinConcentrations() {
    for (Protein p : tfProteins)
      p.concentration = (1.0 - inputConcentration) / tfProteins.length;

    for (Protein p : pProteins)
      p.concentration = 1.0 / pProteins.length;
  }

  /**
   * Normalises the TF protein concentrations so the sum total is 1.0.
   * Input TF proteins are not affected.
   */
  private void normaliseTFProteinConcentrations() {
    double total = 0;
    for (int i = 0; i < tfProteins.length-numberOfInputs; i++)
      total += tfProteins[i].concentration;

    if (total > 0.0)
      for (int i = 0; i < tfProteins.length; i++)
        if (i  < tfProteins.length-numberOfInputs) {
          tfProteins[i].concentration *= 1.0 - inputConcentration;
          tfProteins[i].concentration /= total;
        }
  }

  /**
   * Checks if the model is at rest. For each protein, check if its
   * concentration has changed more than REST_EPSILON over the last
   * REST_STEP time steps.
   *
   * @param results the model data
   * @param t the current timestep
   * @result whether the model has leveled out or not
   */
  private boolean atRest(double[][] results, int t) {
    /* Wait until at least REST_STEP time steps have passed */
    if (t < REST_STEP)
      return false;

    /* Check TF proteins */
    for (int i = 0; i < tfProteins.length; i++)
      if (Math.abs(results[t][i] - results[t-REST_STEP][i]) > REST_EPSILON) {
        return false;
      }

    /* Check P proteins */
    for (int i = 0; i < pProteins.length; i++)
      if (Math.abs(results[t][tfProteins.length+i] -
                   results[t-REST_STEP][tfProteins.length+i]) > REST_EPSILON) {
        return false;
      }

    return true;
  }

  /**
   * Inject input proteins into the model (replacing the current inputs)
   *
   * @param inputProteins Proteins to replace the current input proteins
   */
  public void injectInputs(ArrayList<Protein> inputProteins) {
    /* If the incorrect size size, then resize */
    if (inputProteins.size() != numberOfInputs)
      tfProteins = Arrays.copyOf(tfProteins, tfProteins.length - numberOfInputs+inputProteins.size());

    /* Add the new inputs */
    numberOfInputs = inputProteins.size();
    for (int p = numberOfInputs; p > 0; p--)
      tfProteins[tfProteins.length - p] = inputProteins.get(p - 1);
    calculateInputConcentration(inputProteins);

    /* Normalise the rest of the TF concnetrations */
    normaliseTFProteinConcentrations();

    /* Regenerate the concentration tables */
    generateTables(); //FIX ME: This should just be for these input proteins, not the entire model
  }

  /**
   * Initialise the model by attempting to reach a steady/stable state.
   * Run until at rest or a timestep of INIT_PERIOD has been reached. 
   * The true parameter indicates to check for stability. 
   */
  public void init() {
    initResults = run(INIT_PERIOD, true);
  }

  /**
   * Iterate the network. 
   * 
   * 
   *
   *
   * @param timeSteps the length of the synchronisation step, i.e., how many iterations/timesteps
   * @param initialising whether to check for, and halt on a stable state
   * @return a new array containing the concentration values of all proteins at each timestep.
   */
  public double[][] run(int timeSteps, boolean initialising) {
    /* Results array */
    results = new double[timeSteps+1][tfProteins.length+pProteins.length];

    /* Iterate the network */
    int t;
    for (t = 0; t < timeSteps && (initialising ? !atRest(results, t - 1) : true); t++) {

      /* Record the current state */
      for (int i = 0; i < tfProteins.length; i++)
        results[t][i] = tfProteins[i].concentration;
      for (int i = 0; i < pProteins.length; i++)
        results[t][tfProteins.length+i] = pProteins[i].concentration;

      /* Calculate production rates */
      double[] geneProductionRates = new double[tfGenes.length];
      double[] pGeneProductionRates = new double[pGenes.length];
      for (int i = 0; i < tfGenes.length;  i++) {
        tfGenes[i].index = i;
        geneProductionRates[i] = ProteinProducer.produce(tfGenes[i], tfProteins[i], tfProteins, umax);
      }
      for (int i = 0; i < pGenes.length;  i++) {
        pGenes[i].index = tfGenes.length+i;
        pGeneProductionRates[i] = ProteinProducer.p_produce(pGenes[i], pProteins[i], tfProteins, umax);
      }

      /* Update protein concentrations c += dc/dt */
      for (int i = 0; i < tfGenes.length;  i++) {
        tfProteins[i].concentration += geneProductionRates[i];
        if (tfProteins[i].concentration < ZERO)
          tfProteins[i].concentration = ZERO;
      }
      for (int i = 0; i < pGenes.length;  i++) {
        pProteins[i].concentration += pGeneProductionRates[i];
        if (pProteins[i].concentration < ZERO)
          pProteins[i].concentration = ZERO;
      }

      /* Re-normalise TF concentration levels */
      double total = 0;
      for (int i = 0; i < tfProteins.length-numberOfInputs; i++)
        total += tfProteins[i].concentration;

      if (total > 0.0)
        for (int i = 0; i < tfProteins.length; i++) {
          if (i  < tfProteins.length-numberOfInputs) {
            tfProteins[i].concentration *= 1.0 - inputConcentration;
            tfProteins[i].concentration /= total;
          }
      }

      /* Re-normalise P concentration levels */
      total = 0;
      for (Protein p : pProteins)
        total += p.concentration;

      if (total > 0.0)
        for (int i = 0; i < pProteins.length; i++) {
          pProteins[i].concentration /= total;
        }
    }

    //Record the final state
    for (int i = 0; i < tfProteins.length; i++)
      results[t][i] = tfProteins[i].concentration;
    for (int i = 0; i < pProteins.length; i++)
      results[t][tfProteins.length+i] = pProteins[i].concentration;

    return results;
  }
  
  /**
   * Constructs a random GRN and runs it for a sync of 2000 time steps without inputs.
   * @param args Optionally contains an integer seed
   */
  public static void main(String[] args) {
    int seed = args.length > 0 ? Integer.parseInt(args[0]) : 42;
    int syncSize = args.length > 1 ? Integer.parseInt(args[1]) : 2000;

    Random r = new Random(seed);
    int[] codons = new int[128];
    for (int i = 0; i <  128; i++)
      codons[i] = r.nextInt();

    Grn grn = new Grn(codons);
    grn.pGenes = new Gene[0];
    grn.pProteins = new Protein[0];

    GRNPrinter.printGRNToFile(seed+".grn", grn, grn.run(syncSize, false));
  }

  private static int intFromByteArray(byte[] bytes) {
    return (bytes[0] & 0x000000FF) << 24 |
      (bytes[1] & 0x000000FF) << 16 |
      (bytes[2] & 0x000000FF) << 8 |
      (bytes[3] & 0x000000FF);
  }

  /**
   * Return an integer array of the integer representations of each
   * gene and its components.
   *
   * @return array of codon values encoding all genes
   */
  public int[] getGRNEncoding() {
    int[] codons = new int[(tfGenes.length + pGenes.length) * Gene.SIZE];
    int i = 0;
    for (Gene g : tfGenes) {
      codons[i++] = g.enhancer;
      codons[i++] = g.inhibitor;
      codons[i++] = g.promoter;
      for (int j = 0; j < 5; j++)
        codons[i++] = g.codons[j];
    }

    for (Gene g : pGenes) {
      codons[i++] = g.enhancer;
      codons[i++] = g.inhibitor;
      codons[i++] = g.promoter;
      for (int j = 0; j < 5; j++)
        codons[i++] = g.codons[j];
    }
    return codons;
  }
}
