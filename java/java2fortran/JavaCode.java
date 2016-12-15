// JavaCode.java
class JavaCode {
    final static int MAXSIZE = 10;

    private native int sumsquaredc(int arr[]);

    public static void main(String args[]) {
        System.out.println("-- We are in the Java program JavaCode --");
        JavaCode c = new JavaCode();
        int arr[] = new int[MAXSIZE];
        System.out.println("Initialize the array arr[]");
        for (int i=0; i<MAXSIZE; i++) {
            arr[i] = i;
            System.out.println(i + "   " + arr[i]);
        }
        System.out.println("Call the C code");
        int sum = c.sumsquaredc(arr);
        System.out.println("-- We are back in Java --");
        System.out.println("Contents of arr[]");
        for (int i=0; i<MAXSIZE; i++)
            System.out.println(i + "   " + arr[i]);
        System.out.println("Sum of squares in arr[] = " + sum);
        System.out.println("Exit Java");
    }

    static {
        // Call up the static library libmycodeinc.so created from ccode.c
        System.loadLibrary("mycodeinc");
    }
}
