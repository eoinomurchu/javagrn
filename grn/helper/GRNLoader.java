package grn.helper;

import grn.Gene;
import grn.Grn;
import grn.Protein;
import grn.ProteinProducer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class GRNLoader {

  private GRNLoader() {
    ;
  }

  /**
   * Constructs a Grn object from a genes file.
   *
   * A .genes file consists of a single gene per line. Each line is
   * made up of eight integers, each seperated by a single
   * space. These integers are the 32 bit representation of the gene's
   * enhancer and inhibitor regulatory sites, promotor site and the
   * 160 bit long gene information split into five equal length
   * segments.
   *
   * @param filename The .gene file to open
   * @return the constructed grn
   */
  public static Grn readFromGenesFile(String filename) {
    ArrayList<Gene> tfGenes = new ArrayList<Gene>();
    ArrayList<Gene> pGenes = new ArrayList<Gene>();
    Grn grn = null;

    try {
      FileInputStream fs = new FileInputStream(new File(filename));
      Scanner scanner = new Scanner(fs);

      StringBuffer sb = new StringBuffer();
      if (scanner.hasNext())
        sb.append(scanner.next());
      while (scanner.hasNext())
        sb.append(" "+scanner.next());
      grn = new Grn(sb.toString(), new ArrayList<Protein>());

      fs.close();
    }
    catch (IOException e) {
      e.printStackTrace();
    }

    return grn;
  }

  public static Grn readFromBinaryTextFile(String filename) {
    ArrayList<Integer> genes = new ArrayList<Integer>();
    Grn grn = null;

    try {
      FileInputStream fs = new FileInputStream(new File(filename));
      byte[] chars = new byte[32];
      int len;
      while((len = fs.read(chars)) > 0) {
        genes.add(binaryArrayToInt(chars, len));
      }

      grn = new Grn(convertIntegers(genes));

      fs.close();
    }
    catch (IOException e) {
      e.printStackTrace();
    }

    return grn;
  }

  private static int binaryArrayToInt(byte[] bits, int length) {
    int a = 0;
    for (int i = 0; i < length; i++) {
      if (bits[i] == '0' || bits[i] == '1')
        a |= ((int)bits[i]-48) << (31 - i);
    }
    return a;
  }

  private static int[] convertIntegers(List<Integer> integers)
  {
    int[] ret = new int[integers.size()];
    Iterator<Integer> iterator = integers.iterator();
    for (int i = 0; i < ret.length; i++)
      {
        ret[i] = iterator.next().intValue();
      }
    return ret;
  }
}
