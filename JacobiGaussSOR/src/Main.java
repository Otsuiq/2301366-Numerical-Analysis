import java.util.Scanner;

public class Main {
    public static void main(String[] args) {
        System.out.print("ใส่ Method ที่คุณต้องการ\nJacobi Iterative พิม 1\nGauss-Seidel Iterative Method พิม 2\nSOR พิม 3\nMethod ที่คุณต้องการคือ ");
        int inputMethod = scanner.nextInt();
        System.out.println("-----------------------------------------");
        if(inputMethod == 1){
            System.out.print("ขนาดเมทริกซ์ nxn, n : ");
            int inputDimensionMetrix = scanner.nextInt();
            System.out.print("TOL : ");
            double inputTOL = scanner.nextDouble();
            System.out.print("Loop : ");
            int inputMax = scanner.nextInt();
            System.out.println("-----------------------------------------");
            double[][] metrix = generateMatrix(inputDimensionMetrix, inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double[] b = generateMatrix(inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double[] initialGuess = generateMatrixX0(inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double[] result = jacobiIteration(metrix, b, initialGuess, inputTOL, inputMax);

            System.out.print("Solution: ");
            int i = 1;
            for (double value : result) {
                System.out.print("x" + i +" = " + value + " ");
                i++;
            }
        } else if(inputMethod == 2){
            System.out.print("ขนาดเมทริกซ์ nxn, n : ");
            int inputDimensionMetrix = scanner.nextInt();
            System.out.print("TOL : ");
            double inputTOL = scanner.nextDouble();
            System.out.print("Loop : ");
            int inputMax = scanner.nextInt();
            System.out.println("-----------------------------------------");
            double[][] metrix = generateMatrix(inputDimensionMetrix, inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double[] b = generateMatrix(inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double[] initialGuess = generateMatrixX0(inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double[] result = gaussSeidel(metrix, b, initialGuess, inputTOL, inputMax);

            System.out.print("Solution: ");
            int i = 1;
            for (double value : result) {
                System.out.print("x" + i +" = " + value + " ");
                i++;
            }
        } else if(inputMethod == 3){
            System.out.print("ขนาดเมทริกซ์ nxn ของคุณคือ ");
            int inputDimensionMetrix = scanner.nextInt();
            System.out.print("ใส่ ω ที่คุณต้องการ : ");
            double inputOmega = scanner.nextDouble();
            System.out.print("TOL : ");
            double inputTOL = scanner.nextDouble();
            System.out.print("คุณต้องการวน Loop ไม่เกินกี่รอบ : ");
            int inputMax = scanner.nextInt();
            System.out.println("-----------------------------------------");
            double[][] metrix = generateMatrix(inputDimensionMetrix, inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double[] b = generateMatrix(inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double[] initialGuess = generateMatrixX0(inputDimensionMetrix);
            System.out.println("-----------------------------------------");
            double omega = inputOmega;

            double[] result = sor(metrix, b, initialGuess, omega, inputTOL, inputMax);

            System.out.print("Solution: ");
            int i = 1;
            for (double value : result) {
                System.out.print("x" + i +" = " + value + " ");
                i++;
            }
        }
    }
    public static double[] jacobiIteration(double[][] A, double[] b, double[] initialGuess, double tolerance, int maxIterations) {
        int n = initialGuess.length;
        double[] xNew = new double[n];

        for (int iteration = 0; iteration < maxIterations; iteration++) {
            for (int i = 0; i < n; i++) {
                double sigma = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sigma += A[i][j] * initialGuess[j];
                    }
                }
                xNew[i] = (b[i] - sigma) / A[i][i];
            }

            if (isConverged(xNew, initialGuess, tolerance)) {
                return xNew;
            }

            System.arraycopy(xNew, 0, initialGuess, 0, n);
        }

        throw new RuntimeException("Jacobi Iteration did not converge within the specified number of iterations.");
    }
    public static double[] gaussSeidel(double[][] A, double[] b, double[] initialGuess, double tolerance, int maxIterations) {
        int n = initialGuess.length;
        double[] x = new double[n];
        double[] xNew = new double[n];

        for (int iteration = 0; iteration < maxIterations; iteration++) {
            for (int i = 0; i < n; i++) {
                double sigma1 = 0.0;
                for (int j = 0; j < i; j++) {
                    sigma1 += A[i][j] * xNew[j];
                }

                double sigma2 = 0.0;
                for (int j = i + 1; j < n; j++) {
                    sigma2 += A[i][j] * x[j];
                }

                xNew[i] = (b[i] - sigma1 - sigma2) / A[i][i];
            }

            if (isConverged(xNew, x, tolerance)) {
                return xNew;
            }

            System.arraycopy(xNew, 0, x, 0, n);
        }

        throw new RuntimeException("Gauss-Seidel did not converge within the specified number of iterations.");
    }
    public static double[] sor(double[][] A, double[] b, double[] initialGuess, double omega, double tolerance, int maxIterations) {
        int n = initialGuess.length;
        double[] x = new double[n];
        double[] xNew = new double[n];

        for (int iteration = 0; iteration < maxIterations; iteration++) {
            for (int i = 0; i < n; i++) {
                double sigma1 = 0.0;
                for (int j = 0; j < i; j++) {
                    sigma1 += A[i][j] * xNew[j];
                }

                double sigma2 = 0.0;
                for (int j = i + 1; j < n; j++) {
                    sigma2 += A[i][j] * x[j];
                }

                xNew[i] = (1 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sigma1 - sigma2);
            }

            if (isConverged(xNew, x, tolerance)) {
                return xNew;
            }

            System.arraycopy(xNew, 0, x, 0, n);
        }

        throw new RuntimeException("SOR did not converge within the specified number of iterations.");
    }
    private static double[][] generateMatrix(int rows, int cols) {
        double[][] matrix = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                int row = i + 1 ;
                int col = j + 1 ;
                Scanner scanner = new Scanner(System.in);
                System.out.print("ค่า A ในแถวที่ " + row + " หลัก " + col + " ของคุณคือ ");
                double inputMetrix = scanner.nextDouble();
                matrix[i][j] = inputMetrix;
            }
        }
        return matrix;
    }
    private static double[] generateMatrix(int rows) {
        Scanner scanner = new Scanner(System.in);
        double[] matrix = new double[rows];
        for (int i = 0; i < rows; i++) {
            int row = i + 1 ;
            System.out.print("ค่า b ในแถวที่ " + row + " ของคุณคือ ");
            double inputB = scanner.nextDouble();
            matrix[i] = inputB;
        }
        return matrix;
    }

    private static double[] generateMatrixX0(int rows) {
        Scanner scanner = new Scanner(System.in);
        double[] matrix = new double[rows];
        for (int i = 0; i < rows; i++) {
            int row = i + 1 ;
            System.out.print("ค่า x" + row + " ของคุณคือ ");
            double inputB = scanner.nextDouble();
            matrix[i] = inputB;
        }
        return matrix;
    }
    private static boolean isConverged(double[] xNew, double[] x, double tolerance) {
        for (int i = 0; i < x.length; i++) {
            if (Math.abs(xNew[i] - x[i]) > tolerance) {
                return false;
            }
        }
        return true;
    }
}