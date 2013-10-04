package grn.helper;

import java.util.Arrays;

public class ArrayUtils {

  private ArrayUtils() {
    ;
  }

  /**
   * Fits oversized array to the correct size.
   * Used for gene arrays.
   * 
   * @param array the array to be resized to fit its contents
   * @param n the number of elements in the array
   * @return A copy of the array resize to n
   */
  public static <T> T[] fitArray(T[] array, int n) {
    if (array.length == 0 || n == 0 || array[0] == null)
        return Arrays.copyOf(array, 0);

    return Arrays.copyOf(array, n);
  }

  /** 
   * Checks if the array is full, if so resize to twice the size and
   * return new array.  Otherwise, return the same array
   *
   * @param array to be checked
   * @param n size of array
   * @return The original array, or a copy with double the length
   */
  public static <T>  T[] checkAndResizeArray(T[] array, int n) {
    if (array.length == 0 || n == 0)
      return array;

    return n == array.length ? Arrays.copyOf(array, n*2) : array;
  }
}
