using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace course_project
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            //Метод Крамера -размерность системы от 2 до 3 алгебраических уравнений -определитель не должен быть равен 0 -в консоль должны записываться только численные выражения
            richTextBox1.Text += Environment.NewLine + " ------Метод Крамера------" + Environment.NewLine;
            int n = 3; /* количество уравнений */
            double[,] A = new double[n, n]; /* матрица системы aMain*/
            double[] b = new double[n]; /* вектор правых частей freeVar*/
            double[] x = new double[n]; /* вектор решения answer*/
            //Матрица А
            A[0, 0] = Convert.ToDouble(textBox1.Text);
            A[0, 1] = Convert.ToDouble(textBox2.Text);
            A[0, 2] = Convert.ToDouble(textBox3.Text);
            A[1, 0] = Convert.ToDouble(textBox5.Text);
            A[1, 1] = Convert.ToDouble(textBox6.Text);
            A[1, 2] = Convert.ToDouble(textBox7.Text);
            A[2, 0] = Convert.ToDouble(textBox9.Text);
            A[2, 1] = Convert.ToDouble(textBox10.Text);
            A[2, 2] = Convert.ToDouble(textBox11.Text);
            //Матрица b
            b[0] = Convert.ToDouble(textBox4.Text);
            b[1] = Convert.ToDouble(textBox8.Text);
            b[2] = Convert.ToDouble(textBox12.Text);
            double det(int N, double[,] B)
            {
                //метод вычисляющий определитель матрицы
                if (N == 2) {
                    return B[0, 0] * B[1, 1] - B[0, 1] * B[1, 0];
                }
                return B[0, 0] * (B[1, 1] * B[2, 2] - B[1, 2] * B[2, 1]) - B[0, 1] * (B[1, 0] * B[2, 2] - B[1, 2] * B[2, 0]) +
                B[0, 2] * (B[1, 0] * B[2, 1] - B[1, 1] * B[2, 0]);
            }
            void equal(int N, double[,] a, double[,] B)
            {
                //– метод присваивающий матрицы ( ), где n-размерность матриц
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++) {
                        a[i, j] = B[i, j];
                    }
            }
            void change(int N, int N1, double[,] a, double[] B)
            {
                for (int i = 0; i < N; i++) {
                    a[i, N1] = B[i];
                }

            }
            int SLAU_kramer(int N, double[,] a, double[] B, double[] X)
            {
                //метод реализующий метод Крамера
                double[,] An = new double[3, 3];
                double det1 = det(N, a);
                if (det1 == 0) return 1;
                for (int i = 0; i < N; i++)
                {
                    equal(N, An, a);
                    change(N, i, An, B);
                    X[i] = det(N, An) / det1;
                }
                return 0;
            }
            if (SLAU_kramer(n, A, b, x) == 1)
            {
                richTextBox1.Text += "Система не имеет решение" + Environment.NewLine;
                return;
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    int t = i + 1;
                    richTextBox1.Text += " x" + t + " = " + x[i]+"; " + Environment.NewLine;
                }
            }
        }

        private void button2_Click(object sender, EventArgs e) {
            //Метод Гаусса
            richTextBox1.Text += Environment.NewLine + " ------Метод Гаусса------" + Environment.NewLine;
            int N = Convert.ToInt32(textBoxN.Text);
            double[,] Arr = new double[N, N];
            double[] Arr1 = new double[N];
            double[] X = new double[N];
            //Матрица А
            Arr[0, 0] = Convert.ToDouble(textBox1.Text);
            Arr[0, 1] = Convert.ToDouble(textBox2.Text);
            Arr[0, 2] = Convert.ToDouble(textBox3.Text);
            Arr[1, 0] = Convert.ToDouble(textBox5.Text);
            Arr[1, 1] = Convert.ToDouble(textBox6.Text);
            Arr[1, 2] = Convert.ToDouble(textBox7.Text);
            Arr[2, 0] = Convert.ToDouble(textBox9.Text);
            Arr[2, 1] = Convert.ToDouble(textBox10.Text);
            Arr[2, 2] = Convert.ToDouble(textBox11.Text);
            //Матрица b
            Arr1[0] = Convert.ToDouble(textBox4.Text);
            Arr1[1] = Convert.ToDouble(textBox8.Text);
            Arr1[2] = Convert.ToDouble(textBox12.Text);
            double E = Convert.ToDouble(textBoxE.Text);

            sysout(Arr, Arr1, N);
            X = gauss(Arr, Arr1, N);
            for (int i = 0; i < N; i++)
                richTextBox1.Text += "x[" + (i + 1) + "]=" + X[i] + Environment.NewLine;
            void sysout(double[,] a, double[] y, int n)
            {
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        richTextBox1.Text += a[i, j] + "*x" + j;
                        if (j < n - 1)
                            richTextBox1.Text += " + ";
                    }
                    richTextBox1.Text += " = " + y[i] + Environment.NewLine;
                }
                return;
            }
            double[] gauss(double[,] a, double[] y, int n)
            {
                double[] x = new double[n];
                double max;
                int k, index;
                double eps = Convert.ToDouble(textBoxE.Text); ;  // точность
                k = 0;
                while (k < n)
                {
                    // Поиск строки с максимальным a[i][k]
                    max = Math.Abs(a[k, k]);
                    index = k;
                    for (int i = k + 1; i < n; i++)
                    {
                        if (Math.Abs(a[i, k]) > max)
                        {
                            max = Math.Abs(a[i, k]);
                            index = i;
                        }
                    }
                    // Перестановка строк
                    if (max < eps)
                    {
                        // нет ненулевых диагональных элементов
                        richTextBox1.Text += "Решение получить невозможно из-за нулевого столбца ";
                        richTextBox1.Text += index + " матрицы A";
                        break;
                    }
                    double temp;
                    for (int j = 0; j < n; j++)
                    {
                        temp = a[k, j];
                        a[k, j] = a[index, j];
                        a[index, j] = temp;
                    }
                    temp = y[k];
                    y[k] = y[index];
                    y[index] = temp;
                    // Нормализация уравнений
                    for (int i = k; i < n; i++)
                    {
                        temp = a[i, k];
                        if (Math.Abs(temp) < eps) continue; // для нулевого коэффициента пропустить
                        for (int j = 0; j < n; j++)
                            a[i, j] = a[i, j] / temp;
                        y[i] = y[i] / temp;
                        if (i == k) continue; // уравнение не вычитать само из себя
                        for (int j = 0; j < n; j++)
                            a[i, j] = a[i, j] - a[k, j];
                        y[i] = y[i] - y[k];
                    }
                    k++;
                }
                // обратная подстановка
                for (k = n - 1; k >= 0; k--)
                {
                    x[k] = y[k];
                    for (int i = 0; i < k; i++)
                        y[i] = y[i] - a[i, k] * x[k];
                }
                return x;
            }
            //end
        }

        private void button4_Click(object sender, EventArgs e)
        {
            //Метод простых итераций
            richTextBox1.Text = Environment.NewLine + " ------Метод простой итерации------" + Environment.NewLine;
            int N = Convert.ToInt32(textBoxN.Text);
            double[,] Arr = new double[N, N];
            double[] Arr1 = new double[N];
            //Матрица А
            Arr[0, 0] = Convert.ToDouble(textBox1.Text);
            Arr[0, 1] = Convert.ToDouble(textBox2.Text);
            Arr[0, 2] = Convert.ToDouble(textBox3.Text);
            Arr[1, 0] = Convert.ToDouble(textBox5.Text);
            Arr[1, 1] = Convert.ToDouble(textBox6.Text);
            Arr[1, 2] = Convert.ToDouble(textBox7.Text);
            Arr[2, 0] = Convert.ToDouble(textBox9.Text);
            Arr[2, 1] = Convert.ToDouble(textBox10.Text);
            Arr[2, 2] = Convert.ToDouble(textBox11.Text);
            //Матрица b
            Arr1[0] = Convert.ToDouble(textBox4.Text);
            Arr1[1] = Convert.ToDouble(textBox8.Text);
            Arr1[2] = Convert.ToDouble(textBox12.Text);
            double E = Convert.ToDouble(textBoxE.Text);
            double[] arr1out;
            Show(Arr);
            Arr = Sys(Arr, Arr1, out arr1out);
            richTextBox1.Text += Environment.NewLine;
            Show(Arr);
            richTextBox1.Text += Environment.NewLine;
            Show1(arr1out);
            richTextBox1.Text += Environment.NewLine;
            Approximation(Arr, arr1out, E, 1);
            richTextBox1.Text += Environment.NewLine;

            double[,] Sys(double[,] matrix, double[] arr1, out double[] arr2) {
                double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)];
                arr2 = new double[arr1.GetLength(0)];
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    arr2[i] = arr1[i] / matrix[i, i];
                    for (int j = 0; j < matrix.GetLength(1); j++)
                    {
                        if (i != j)
                            result[i, j] = -(matrix[i, j] / matrix[i, i]);
                        else
                            result[i, j] = 0;
                    }
                }
                return result;
            }
            void Approximation(double[,] matrix, double[] b, double Et, int n)
            {
                int count = 0;
                double[] x = new double[b.GetLength(0)];
                double[] x0 = new double[b.GetLength(0)];
                double sum = 0;
                for (int root = 0; root < b.GetLength(0); root++)
                {
                    count = 0;
                    for (int i = 0; i < b.GetLength(0); i++)
                    {
                        x0[i] = b[i];
                    }

                    for (int i = 0; i < matrix.GetLength(0); i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            sum = matrix[i, j] * x0[j] + sum;
                        }
                        x[i] = b[i] + sum;
                        sum = 0;
                    }
                    while (Math.Abs(x[root] - x0[root]) >= Et)
                    {
                        for (int k = 0; k < matrix.GetLength(0); k++)
                        {
                            x0[k] = x[k];
                        }
                        for (int i = 0; i < matrix.GetLength(0); i++)
                        {
                            for (int j = 0; j < 3; j++)
                            {
                                sum = matrix[i, j] * x0[j] + sum;
                            }
                            x[i] = b[i] + sum;
                            sum = 0;
                        }
                        richTextBox1.Text += Environment.NewLine + x[root];
                        count++;
                    }
                    richTextBox1.Text += Environment.NewLine + "x[" + root + "]=";
                    richTextBox1.Text += Environment.NewLine + x[root];
                    richTextBox1.Text += Environment.NewLine + count;
                }
            }
            void Show(double[,] matrix)
            {
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    for (int j = 0; j < matrix.GetLength(1); j++)
                    {
                        richTextBox1.Text += "\t" + matrix[i, j] + "\t";
                    }
                    richTextBox1.Text += Environment.NewLine;
                }
            }
            void Show1(double[] matrix)
            {
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    richTextBox1.Text += "\t";
                    int t = i + 1;
                    richTextBox1.Text += Environment.NewLine +"x" + t + " = " + matrix[i];
                }
            }

            //end
        }

        private void button6_Click(object sender, EventArgs e)
        {
            //метод Бройдена
            richTextBox1.Text += Environment.NewLine + " ------Метод Бройдена------" + Environment.NewLine;
            int N = Convert.ToInt32(textBoxN.Text); ;
            double[,] yakob = new double[N, N];
            yakob[0, 0] = Convert.ToDouble(textBox1.Text);
            yakob[0, 1] = Convert.ToDouble(textBox2.Text);
            yakob[0, 2] = Convert.ToDouble(textBox3.Text);
            yakob[1, 0] = Convert.ToDouble(textBox5.Text);
            yakob[1, 1] = Convert.ToDouble(textBox6.Text);
            yakob[1, 2] = Convert.ToDouble(textBox7.Text);
            yakob[2, 0] = Convert.ToDouble(textBox9.Text);
            yakob[2, 1] = Convert.ToDouble(textBox10.Text);
            yakob[2, 2] = Convert.ToDouble(textBox11.Text);
            double[] V = new double[N];
            double[] B = new double[N];
            double[] Bnach = new double[N];
            double E = Convert.ToDouble(textBoxE.Text); //удовлетворяющую погрешность
            int maunI = 0;
            int naid = 0;
            int stop = 0;
            double S=0;
            richTextBox1.Text += "Матрица Якоби:" + Environment.NewLine;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    richTextBox1.Text += yakob[i, j] + " ";
                }
                richTextBox1.Text += Environment.NewLine;
            }
            while ((maunI != 10) && (naid != 1) && (stop != 1))
            {
                maunI++;
                Bnach[0] = V[0] + V[1] - 3;
                Bnach[1] = V[0] * V[0] + V[1] * V[1] - 9;
                int iter = 0;
                double[,] A = new double[N, N];
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        A[i, j] = yakob[i, j];
                    }
                }
                while (iter != N - 1)
                {
                    for (int h = 0; h < N; h++) { B[h] = Bnach[h] * (-1); }
                    double pomny = A[iter, iter];
                    for (int j = iter; j < N; j++)
                    {
                        A[iter, j] = A[iter, j] / pomny;
                    }
                    B[iter] = B[iter] / pomny;
                    for (int i = iter + 1; i < N; i++)
                    {
                        double zap = A[i, iter];
                        for (int j = iter; j < N; j++)
                        {
                            A[i, j] = A[i, j] - A[iter, j] * zap;
                        }
                        B[i] = B[i] - B[iter] * zap;
                    }
                    iter++;
                }
                double[] X = new double[N];
                if (A[N - 1, N - 1] != 0)
                {
                    X[N - 1] = B[N - 1] / A[N - 1, N - 1];
                }
                else X[N - 1] = 0;
                double SYM = 0;
                int l = N - 2;
                for (int i = N - 2; i >= 0; i--)
                {
                    SYM = 0;
                    for (int j = i + 1; j <= N - 1; j++)
                    {
                        SYM = SYM + A[i, j] * X[j];
                    }
                    if (A[i, l] != 0)
                    {
                        X[i] = (B[i] - SYM) / A[i, l];
                    }
                    else X[i] = 0;
                    l--;
                }
                double[] XJ = new double[N];
                double promq = 0; double mq = 0; double nq = 0; S = 0;
                for (int i = 0; i < N; i++)
                {
                    XJ[i] = V[i] + X[i];
                    if (X[i] >= 0)
                    {
                        promq = X[i] + promq;
                    }
                    else
                    {
                        promq = -X[i] + promq;
                    }
                    if (V[i] >= 0)
                    {
                        mq = mq + V[i];
                    }
                    else
                    {
                        mq = mq - V[i];
                    }
                    if (XJ[i] >= 0)
                    {
                        nq = nq + XJ[i];
                    }
                    else
                    {
                        nq = nq - XJ[i];
                    }
                }
                if (mq != 0)
                {
                    S = promq / mq;
                }
                else
                {
                    S = promq / nq;
                }
                if (S < 0)
                {
                    S = -S;
                }
                if (S < E)
                {
                    richTextBox1.Text += "S " + S + " " + Environment.NewLine;
                    naid = 1;
                    richTextBox1.Text += "Найдено решение: "+ Environment.NewLine;
                    for (int i = 0; i < N; i++)
                    {
                        richTextBox1.Text += XJ[i] + " " + Environment.NewLine;
                    }
                    richTextBox1.Text += " Количество итераций " + maunI;
                }
                else
                {
                    if (S > 20)
                    {
                        richTextBox1.Text += " Процесс расходится "; stop = 1;
                    }
                    else
                    {
                        if (maunI == 10)
                        {
                            richTextBox1.Text += " За 10 итераций решение не найдено ";
                        }
                        else
                        {
                            double[] Y = new double[N];
                            Y[0] = (XJ[0] + XJ[1] - 3) - Bnach[0];
                            Y[1] = (XJ[0] * XJ[0] + XJ[1] * XJ[1] - 9) - Bnach[1];
                            double[,] J = new double[N, N];
                            for (int i = 0; i < N; i++)
                            {
                                for (int j = 0; j < N; j++)
                                {
                                    J[i, j] = yakob[i, j];
                                    yakob[i, j] = 0;
                                }
                            }
                            double[] ymnMAS = new double[N]; double[] PRMAS = new double[N];
                            for (int i = 0; i < N; i++)
                            {
                                double Ymn = 0;
                                for (int j = 0; j < N; j++)
                                {
                                    Ymn = Ymn + J[i, j] * X[j];
                                }
                                ymnMAS[i] = Ymn;
                                PRMAS[i] = Y[i] - ymnMAS[i];
                            }
                            double del = 0;
                            for (int i = 0; i < N; i++)
                            {
                                del = del + X[i] * X[i];
                            }
                            for (int i = 0; i < N; i++)
                            {
                                for (int j = 0; j < N; j++)
                                {
                                    yakob[i, j] = J[i, j] + ((PRMAS[i] * X[j]) / del);
                                }
                            }
                            for (int i = 0; i < N; i++)
                            {
                                V[i] = XJ[i];
                            }
                        }
                    }
                }
            }
                        
        }

        private void button7_Click(object sender, EventArgs e)
        {
            //Метод вращений
            richTextBox1.Text += Environment.NewLine + " ------Метод вращений------" + Environment.NewLine;
            int N = Convert.ToInt32(textBoxN.Text) + 1;
            double[,] a = new double[N, N]; //a1-промежуточный массив
            double[] B = new double[N];
            //Матрица А
            a[0, 0] = Convert.ToDouble(textBox1.Text);
            a[0, 1] = Convert.ToDouble(textBox2.Text);
            a[0, 2] = Convert.ToDouble(textBox3.Text);
            a[1, 0] = Convert.ToDouble(textBox5.Text);
            a[1, 1] = Convert.ToDouble(textBox6.Text);
            a[1, 2] = Convert.ToDouble(textBox7.Text);
            a[2, 0] = Convert.ToDouble(textBox9.Text);
            a[2, 1] = Convert.ToDouble(textBox10.Text);
            a[2, 2] = Convert.ToDouble(textBox11.Text);
            //Матрица B
            B[0] = Convert.ToDouble(textBox4.Text);
            B[1] = Convert.ToDouble(textBox8.Text);
            B[2] = Convert.ToDouble(textBox12.Text);
            Rotations(a, B, N);
            void Rotations(double[,] A, double[] X, int i)
            {
                {
                    double C, S, A1;

                    for (int n = 0; n < i; n++)
                    {
                        C = A[n, n] / (Math.Sqrt((A[n, n] * A[n, n]) + (A[n, n] * A[n, n])));
                        S = A[n, n] / (Math.Sqrt((A[n, n] * A[n, n]) + (A[n, n] * A[n, n])));
                        for (int m = 0; m < i; m++)
                        {
                            A1 = A[n, m];
                            A[n, m] = C * A1 + S * A[n, m];
                            A[n, m] = -S * A1 + C * A[n, m];
                        }
                    }
                    for (int n = i - 2; n >= 0; n--)
                    {
                        for (int m = n + 1; m < i; m++)
                        {
                            X[n] -= A[n, m] * X[m];
                        }
                        X[n] /= A[n, n];
                        richTextBox1.Text += Environment.NewLine + " x " + (n + 1) + " = " + X[n];
                    }
                }
            }
        }

        private void button8_Click(object sender, EventArgs e)
        {
            //Метод Зейделя
            richTextBox1.Text = Environment.NewLine + " ------Метод Зейделя------" + Environment.NewLine;
            int N = Convert.ToInt32(textBoxN.Text); //переменная N 
            double E = Convert.ToDouble(textBoxE.Text); //переменная E точности.
            double[,] a = new double[N, N];
            double[] b = new double[N];
            //Матрица a
            a[0, 0] = Convert.ToDouble(textBox1.Text);
            a[0, 1] = Convert.ToDouble(textBox2.Text);
            a[0, 2] = Convert.ToDouble(textBox3.Text);
            a[1, 0] = Convert.ToDouble(textBox5.Text);
            a[1, 1] = Convert.ToDouble(textBox6.Text);
            a[1, 2] = Convert.ToDouble(textBox7.Text);
            a[2, 0] = Convert.ToDouble(textBox9.Text);
            a[2, 1] = Convert.ToDouble(textBox10.Text);
            a[2, 2] = Convert.ToDouble(textBox11.Text);
            //Матрица b
            b[0] = Convert.ToDouble(textBox4.Text);
            b[1] = Convert.ToDouble(textBox8.Text);
            b[2] = Convert.ToDouble(textBox12.Text);
            iterat(a, b, N, E);
            double firstNorm(double[,] A,  int n, int m) {
                int i, j;
                double sum = 0, subSum;
                for (i = 0; i < n; i++)
                {
                    subSum = 0;
                    for (j = 0; j < m; j++)
                    {
                        subSum += Math.Abs(A[i, j]);
                    }
                    if (subSum > sum)
                    {
                        sum = subSum;
                    }
                }
                return sum;
            }
            double secondNorm(double[,] A, int n, int m) {
                int i, j;
                double sum = 0, subSum;
                for (j = 0; j < n; j++)
                {
                    subSum = 0;
                    for (i = 0; i < m; i++)
                    {
                        subSum += Math.Abs(A[i, j]);
                    }
                    if (subSum > sum)
                    {
                        sum = subSum;
                    }
                }
                return sum;
            }
            double thirdNorm(double[,] A, int n, int m) {
                int i, j;
                double sum = 0;
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < m; j++)
                    {
                        sum += (A[i, j] * A[i, j]);
                    }
                }
                sum = Math.Sqrt(sum);
                return sum;
            }

            double okr(double X, double eps)
            {
                int i = 0;
                while (eps != 1)
                {
                    i++;
                    eps *= 10;
                }
                double Okr = Math.Pow(10, i);
                X = (X * Okr + 0.5) / (Okr);
                return X;
            }
            void iterat(double[,] A, double[] B, int n, double eps) {
                int k = 0;
                int i, j;
                double[] X = new double[n];
                double[] Y = new double[n];
                double s;
                double g;

                for (i = 0; i < n; i++)
                {
                    X[i] = B[i];
                }
                do
                {
                    s = 0; k++;
                    { // Решаем систему методом Зейделя.
                        for (i = 0; i < n; i++)
                        {
                            g = B[i];
                            for (j = 0; j < n; j++)
                            {
                                g = g + A[i, j] * X[j];
                            }
                            s += (X[i] - g) * (X[i] - g);
                            X[i] = g;
                            richTextBox1.Text += Environment.NewLine + "x" + (i + 1) + " = " + g;
                        }
                    }
                } while (Math.Sqrt(s) >= eps * (1 - thirdNorm(A, n, n)) / thirdNorm(A, n, n));
                richTextBox1.Text += Environment.NewLine + "Решение системы:";
                for (i = 0; i < n; i++)
                {
                    richTextBox1.Text += Environment.NewLine + "X" + i + " = " + okr(X[i], eps);
                }
                richTextBox1.Text += Environment.NewLine + "Число итераций: " + (k - 1);
                richTextBox1.Text += Environment.NewLine + "Первая норма матрицы A: " + firstNorm(A, n, n);
                richTextBox1.Text += Environment.NewLine + "Вторая норма матрицы A: " + secondNorm(A, n, n);
                richTextBox1.Text += Environment.NewLine + "Третья норма матрицы A: " + thirdNorm(A, n, n);

            }
        }

        private void button9_Click(object sender, EventArgs e)
        {
            //Метод Жордана-Гаусса
            richTextBox1.Text += Environment.NewLine + "Метод Жордана-Гаусса: " + Environment.NewLine;
            int n = Convert.ToInt32(textBoxN.Text);
            int m = Convert.ToInt32(textBoxN.Text);
            m += 1;
            int i, j;
            double[,] v1 = new double[n, m];
            //Матрица a
            v1[0, 0] = Convert.ToDouble(textBox1.Text);
            v1[0, 1] = Convert.ToDouble(textBox2.Text);
            v1[0, 2] = Convert.ToDouble(textBox3.Text);
            v1[0, 3] = Convert.ToDouble(textBox4.Text);
            v1[1, 0] = Convert.ToDouble(textBox5.Text);
            v1[1, 1] = Convert.ToDouble(textBox6.Text);
            v1[1, 2] = Convert.ToDouble(textBox7.Text);
            v1[1, 3] = Convert.ToDouble(textBox8.Text);
            v1[2, 0] = Convert.ToDouble(textBox9.Text);
            v1[2, 1] = Convert.ToDouble(textBox10.Text);
            v1[2, 2] = Convert.ToDouble(textBox11.Text);
            v1[2, 3] = Convert.ToDouble(textBox12.Text);
            richTextBox1.Text += Gauss(n, m, v1);

            string Gauss(int Rows, int Column, double[,] matr)
            {
                int k, q;
                double v;
                string answer = "";
                for (q = 0; q < Rows; q++)
                {
                    //делаем главную диогональ единицами  
                    v = matr[q, q];
                    for (k = 0; k < Column; k++)
                        matr[q, k] /= v;
                    //обнуляем числа под единицами главной диогoнали
                    for (int i1 = q + 1; i1 < Rows; i1++)
                    {
                        v = matr[i1, q];
                        for (k = q; k < Column; k++)
                            matr[i1, k] = matr[i1, k] - matr[q, k] * v;
                    }
                }
                for (q = 0; q < Rows; q++)
                    for (int i2 = 0; i2 < (Rows - 1) - q; i2++)
                    {
                        v = matr[i2, (Column - 1) - q - 1];
                        for (k = Column - 1 - q - 1; k < Column; k++)
                            matr[i2, k] = matr[i2, k] - matr[(Rows - 1) - q, k] * v;
                    }
                for (int i3 = 0; i3 < Rows; i3++)
                    answer += "x_" + i3 + " = " + matr[i3, Column - 1] + "\r\n";
                return answer;
            }
            //end
        }


        private void Form1_Load(object sender, EventArgs e)
        {

        }
    }
}
