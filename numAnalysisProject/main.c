#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printMatrix(int N, float **mat) {
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      printf("%f\t", mat[i][j]);
    }
    printf("\n");
  }
};

typedef struct {
  char **tokens;
  int size;
  int capacity;
} TokenList;

void addToken(TokenList *list, const char *token) {
  if (list->size >= list->capacity) {
    list->capacity *= 2;
    list->tokens = realloc(list->tokens, list->capacity * sizeof(char *));
  }
  list->tokens[list->size++] = strdup(token);
}

TokenList *tokenize(const char *expression) {
  TokenList *tokens = malloc(sizeof(TokenList));
  tokens->size = 0;
  tokens->capacity = 10;
  tokens->tokens = malloc(tokens->capacity * sizeof(char *));

  const char *ptr = expression;
  while (*ptr) {
    if (isdigit(*ptr)) {
      const char *start = ptr;
      while (isdigit(*ptr) || *ptr == '.')
        ptr++;
      int len = ptr - start;
      char *number = malloc(len + 1);
      strncpy(number, start, len);
      number[len] = '\0';
      addToken(tokens, number);
      free(number);
    } else if (strchr("+-*/^()", *ptr)) {
      char op[2] = {*ptr, '\0'};
      addToken(tokens, op);
      ptr++;
    } else if (isalpha(*ptr)) {
      const char *start = ptr;
      while (isalpha(*ptr))
        ptr++;
      int len = ptr - start;
      char *func = malloc(len + 1);
      strncpy(func, start, len);
      func[len] = '\0';
      if (strcmp(func, "log") == 0) {
        addToken(tokens, strdup("l"));
      } else {
          addToken(tokens, func);
      }
      free(func);
    } else {
      ptr++;
    }
  }

  return tokens;
}

void freeTokenList(TokenList *tokens) {
  for (int i = 0; i < tokens->size; i++) {
    free(tokens->tokens[i]);
  }
  free(tokens->tokens);
  free(tokens);
}

typedef struct {
  char **data;
  int top;
  int capacity;
} Stack;

Stack *create_stack(int capacity) {
  Stack *stack = malloc(sizeof(Stack));
  stack->capacity = capacity;
  stack->top = -1;
  stack->data = malloc(stack->capacity * sizeof(char *));
  return stack;
}

void push(Stack *stack, char *item) { // stack e eleman ekle
  stack->data[++stack->top] = strdup(item);
}

char *pop(Stack *stack) { // eleman çıkar
  if (stack->top == -1)
    return NULL;
  char *item = stack->data[stack->top];
  stack->data[stack->top--] = NULL;
  return item;
}

char *peek(Stack *stack) { // çıkarmadan değerini döndürür
  if (stack->top == -1)
    return NULL;
  return stack->data[stack->top];
}

void free_stack(Stack *stack) {
  for (int i = 0; i <= stack->top; i++) {
    free(stack->data[i]);
  }
  free(stack->data);
  free(stack);
}

int precedence(char *op) { // işlem önceliği, öncelikli olan >
  if (strcmp(op, "+") == 0 || strcmp(op, "-") == 0)
    return 1;
  if (strcmp(op, "*") == 0 || strcmp(op, "/") == 0)
    return 2;
  if (strcmp(op, "^") == 0 || strcmp(op, "l") == 0)
    return 3;
  return 0;
}

int is_function(char *token) { // fonksiyon mu
  return strcmp(token, "sin") == 0 || strcmp(token, "cos") == 0 ||
         strcmp(token, "tan") == 0 || strcmp(token, "cot") == 0 ||
         strcmp(token, "arcsin") == 0 || strcmp(token, "arccos") == 0 ||
         strcmp(token, "arctan") == 0 || strcmp(token, "arccot") == 0;
}

int is_variable(char *token) { return strcmp(token, "x") == 0; }

TokenList *infixToPostfix(TokenList *tokens) {
  Stack *stack = create_stack(tokens->size);
  TokenList *output = malloc(sizeof(TokenList));
  output->size = 0;
  output->capacity = tokens->size;
  output->tokens = malloc(output->capacity * sizeof(char *));

  for (int i = 0; i < tokens->size; i++) {
    char *token = tokens->tokens[i];
    if (isdigit(token[0]) || is_variable(token)) {
      output->tokens[output->size++] = strdup(token);
    } else if (strchr("+-*/^l", token[0])) {
      while (stack->top != -1 && precedence(peek(stack)) >= precedence(token)) {
        output->tokens[output->size++] = pop(stack);
      }
      push(stack, token);
    } else if (strcmp(token, "(") == 0) {
      push(stack, token);
    } else if (strcmp(token, ")") == 0) {
      while (stack->top != -1 && strcmp(peek(stack), "(") != 0) {
        output->tokens[output->size++] = pop(stack);
      }
      pop(stack);
      if (stack->top != -1 && is_function(peek(stack))) {
        output->tokens[output->size++] = pop(stack);
      }
    } else if (is_function(token)) {
      push(stack, token);
    }
  }

  while (stack->top != -1) {
    output->tokens[output->size++] = pop(stack);
  }

  free_stack(stack);
  return output;
}

double evaluateFunction(char *func, double value) {
  if (strcmp(func, "sin") == 0)
    return sin(value);
  if (strcmp(func, "cos") == 0)
    return cos(value);
  if (strcmp(func, "tan") == 0)
    return tan(value);
  if (strcmp(func, "cot") == 0)
    return tan(1 / value);
  if (strcmp(func, "arcsin") == 0)
    return asin(value);
  if (strcmp(func, "arccos") == 0)
    return acos(value);
  if (strcmp(func, "arctan") == 0)
    return atan(value);
  if (strcmp(func, "arccot") == 0)
    return atan(1 / value);
  return 0.0;
}

double ValuePostfix(TokenList *postfix, double x) {
  Stack *stack = create_stack(postfix->size);

  for (int i = 0; i < postfix->size; i++) {
    char *token = postfix->tokens[i];
    if (isdigit(token[0]) || (token[0] == '.' && isdigit(token[1]))) {
      push(stack, token);
    } else if (is_variable(token)) {
      double value = x;
      char buffer[50];
      sprintf(buffer, "%lf", value);
      push(stack, buffer);
    } else if (strchr("+-*/^l", token[0])) {
      double b = strtod(pop(stack), NULL);
      double a = strtod(pop(stack), NULL);
      double result;
      if (strcmp(token, "+") == 0)
        result = a + b;
      else if (strcmp(token, "-") == 0)
        result = a - b;
      else if (strcmp(token, "*") == 0)
        result = a * b;
      else if (strcmp(token, "/") == 0)
        result = a / b;
      else if (strcmp(token, "^") == 0)
        result = pow(a, b);
      else if (strcmp(token, "l") == 0)
        result = log(b) / log(a);
      char buffer[50];
      sprintf(buffer, "%lf", result);
      push(stack, buffer);
    } else if (is_function(token)) {
      double value = strtod(pop(stack), NULL);
      double result = evaluateFunction(token, value);
      char buffer[50];
      sprintf(buffer, "%lf", result);
      push(stack, buffer);
    }
  }

  double result = strtod(pop(stack), NULL);
  free_stack(stack);
  return result;
}

double calculateFunction(char *function, double x) {
  TokenList *tokens = tokenize(function); // fonksiyon stringini parçalar
  TokenList *postfix =
      infixToPostfix(tokens); // operand ve operatörlerin postfix sırlamasına(postfix bilgisayarın okuyabilmesi için)
                                // döndürür (operatörler sonda)

  double result = ValuePostfix(postfix, x); // yeni sıralamaya göre hesap

  freeTokenList(tokens);
  freeTokenList(postfix);

  return result;
}

double Bisection() {
  double lower, upper, middle, tolerans;
  double step = 1;

  char functionString[128];

  printf("Kökü bulunacacak denklemi giriniz:");
  scanf("%s", functionString);
  printf("Taranacak aralığı giriniz:");
  scanf("%lf %lf", &lower, &upper);
  printf("Tolerans değerini giriniz:");
  scanf("%lf", &tolerans);

  middle = (lower + upper) / 2;
  while (fabs(lower - upper) / step > tolerans) {
    double u = calculateFunction(functionString, upper);
    double m = calculateFunction(functionString, middle);
    if (u * m < 0) {
      lower = middle;
      middle = (upper + middle) / 2;
    } else {
      upper = middle;
      middle = (lower + middle) / 2;
    }

    step *= 2;
  }

  return middle;
}

double RegulaFalsi() {
  double a, b, tolerans;
  double step = 1;
  char functionString[128];
  int n;

  printf("Kökü bulunacacak denklemi giriniz:");
  scanf("%s", functionString);
  printf("Taranacak aralığı giriniz:");
  scanf("%lf %lf", &a, &b);
  printf("Tolerans değerini giriniz:");
  scanf("%lf", &tolerans);
  if (calculateFunction(functionString, a) *
          calculateFunction(functionString, b) >=
      0) {
    printf("Belirtilen aralıkta kök yok.\n");
    return 0;
  }
  double result = a;
  while (fabs(a - b) / step >= tolerans) {
    double f_a = calculateFunction(functionString, a);
    double f_b = calculateFunction(functionString, b);

    result = (b * f_a - a * f_b) / (f_a - f_b);
    if (calculateFunction(functionString, result) * f_a < 0) {
      b = result;
    } else {
      a = result;
    }

    step *= 2;
  }

  return result;
}

double NewtonRaphson() {
  double tolerans;
  int i = 0;
  double x[3];

  char functionString[128], derivativeString[128];

  printf("Kökü bulunacacak denklemi girin:");
  scanf("%s", functionString);
  printf("Kökü bulunacak denklemin türevini girin:");
  scanf("%s", derivativeString);
  printf("Başlangıç değerini giriniz:");
  scanf("%lf", &x[0]);
  printf("Tolerans değerini giriniz:");
  scanf("%lf", &tolerans);

  do {
    x[(i + 1) % 3] = x[i] - calculateFunction(functionString, x[i]) /
                                calculateFunction(derivativeString, x[i]);
    i = (i + 1) % 3;
  } while (fabs(x[i] - x[(i - 1) % 3]) >= tolerans);

  return x[i];
}

float **createMatrix(int n) {
  int i;
  float **matrix = (float **)malloc(sizeof(float *) * n);
  for (i = 0; i < n; i++) {
    matrix[i] = (float *)malloc(sizeof(float) * n);
  }
  return matrix;
}

void freeMatrix(float **matrix, int n) {
  int i;
  for (i = 0; i < n; i++) {
    free(matrix[i]);
  }
  free(matrix);
}

void inverseMatrix() { // NxN lik matrisin tersi
  int N, i, j, k;
  float eliminationCoefficient, pivotInverse;

  printf("Kare matrisin satır ve sütun sayısını giriniz:");
  scanf("%d", &N);
  float **mat = createMatrix(N);
  float **inverse = createMatrix(N);
  printf("Matris elemanlarini girin:\n");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("[%d][%d]: ", i + 1, j + 1);
      scanf("%f", &mat[i][j]);

      if (i == j) {
        inverse[i][j] = 1.0;
      } else {
        inverse[i][j] = 0.0;
      }
    }
  }
  printf("Matrisin tersi:\n");
  // Matrisin tersini bulmak için Gauss-Jordan yöntemi
  for (i = 0; i < N; i++) {
    if (mat[i][i] == 0.0) { // diyagonalde 0 varsa hesaplanamaz
      printf("Ters matris bulunamadı\n");
      return;
    }
    pivotInverse = 1.0 / mat[i][i];
    for (j = 0; j < N; j++) {
      mat[i][j] *= pivotInverse;
      inverse[i][j] *= pivotInverse;
    }
    for (k = 0; k < N; k++) {
      if (k != i) {
        eliminationCoefficient = mat[k][i];
        for (j = 0; j < N; j++) {
          mat[k][j] -= eliminationCoefficient * mat[i][j];
          inverse[k][j] -= eliminationCoefficient * inverse[i][j];
        }
      }
    }
  }
  printMatrix(N, (float **)inverse);
  freeMatrix(mat, N);
  freeMatrix(inverse, N);
}

void GaussElemination() { // matrisi üst üçgensel forma getirerek çözümü
                          // kolaylaştırır
  float **matrix;
  float scalingFactor; // satırları birbirinden çıkarmak için kullanılır
  int N, i, j, k;
  printf("Kare matrisin satır ve sütun sayısını giriniz:");
  scanf("%d", &N);
  matrix = createMatrix(N + 1);
  printf("Genişletilmiş kat sayılar matrisi girin:\n"); // denklem sistemi
                                                        // genişletişlmiş kat
                                                        // sayılar matrisi
                                                        // şeklinde yazılır
  for (i = 0; i < N; i++) {
    for (j = 0; j < N + 1; j++) {
      printf("Eleman [%d][%d]: ", i + 1, j + 1);
      scanf("%f", &matrix[i][j]);
    }
  }
  for (i = 0; i < N - 1; i++) {
    if (matrix[i][i] == 0.0) {
      printf("Matrsin tersi alınamıyor\n"); // pivot eleman 0 ise matrisin tersi
                                            // alınamaz
      return;
    }
    for (k = i + 1; k < N; k++) {
      scalingFactor = matrix[k][i] / matrix[i][i];
      for (j = 0; j <= N; j++) {
        matrix[k][j] -= scalingFactor * matrix[i][j];
      }
    }
  }
  printf("Üst üçgensel forma dönüştürülmüş matris:\n");
  printMatrix(N, matrix);

  freeMatrix(matrix, N + 1);
}

void GaussSeidal() {
  int N, max, i, j, iteration;
  float Matrix[128][128], MatrixResults[128], x[128], x2[128];

  printf("Matrisin satir ve sütun sayısını girin: ");
  scanf("%d", &N);

  printf("İterasyon sayisini girin: "); // Durma koşulu olarak iterasyon sınırı
                                        // koydum
  scanf("%d", &max);

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      printf("Matrisin [%d][%d]. elemanini girin: ", i + 1, j + 1);
      scanf("%f", &Matrix[i][j]);
    }
  }

  printf("Sonuçları giriniz:\n");
  for (i = 0; i < N; i++) {
    scanf("%f", &MatrixResults[i]);
  }
  printf("Başlangıç değerlerini girin:\n");
  for (i = 0; i < N; i++) {
    printf("x%d: ", i + 1);
    scanf("%f", &x[i]);
  }
  for (iteration = 0; iteration < max; iteration++) {
    for (i = 0; i < N; i++) {
      x2[i] = MatrixResults[i];
      for (j = 0; j < N; j++) {
        if (j != i) {
          x2[i] -= Matrix[i][j] * x[j];
        }
      }
      x2[i] /= Matrix[i][i];
      x[i] = x2[i]; // yeni x değeri günceller. bu yöntemde jacobi iterasyon
                    // metodundan farklı olarak her iterasyonda yeni bulunan
                    // değerler kullanılır.
    }

    printf("%d. İterasyon:\n", iteration + 1);
    for (i = 0; i < N; i++) {
      printf("x%d = %f\n", i + 1, x[i]);
    }
    printf("\n");
  }
}

double sayisalTurev() {
  double x, h, result;
  char funcitonString[128];
  int option;
  printf("Sayısal türevi hesaplanacak fonksiyonu giriniz:");
  scanf("%s", funcitonString);
  printf("İleri farklar yöntemi için: 1\nGeri farklar yöntemi için: 2\nMerkezi "
         "farklar yöntemi için: 3\n");
  scanf("%d", &option);
  printf("x değerini girin: \n");
  scanf("%lf", &x);
  printf("h değerini girin: \n");
  scanf("%lf", &h);

  if (option == 1) { // ileri farklar yöntemi ile türev
    result = (calculateFunction(funcitonString, x + h) -
              calculateFunction(funcitonString, x)) /
             h;

  } else if (option == 2) { // geri farklar yöntemi ile türev
    result = (calculateFunction(funcitonString, x) -
              calculateFunction(funcitonString, x - h)) /
             h;

  } else if (option == 3) { // merkezi farklar yöntemi ile türev
    result = (calculateFunction(funcitonString, x + h) -
              calculateFunction(funcitonString, x - h)) /
             (2 * h);

  } else
    return 0;
  return result;
}

double Simpson() {
  double x, y, result = 0;
  int n, a, b, i, j = 1, option;
  char functionString[128];

  printf("İntegrali hesaplanacak fonksiyonu giriniz:");
  scanf("%s", functionString);
  printf("integralin başlangıç değerini girin: ");
  scanf("%d", &a);
  printf("integralin son değerini girin: ");
  scanf("%d", &b);
  printf("n sayısını girin: "); // integral kaça bölünecek
  scanf("%d", &n);
  printf("Simpson 1/3 icin : 1\nSimpson 3/8 icin : 2\n");
  scanf("%d", &option);
  x = (double)(b - a) / (double)n;
  if (x < 0) {
    printf("Sınırları kontrol edniz");
  } else if (x == 0) {
    result = 0;             // sınır değerler aynı ise işlem yapm gerk yok
  } else if (option == 1) { // simpson 1/3 için
    for (i = 0; i < n; i++) {
      result += (x *
                 (calculateFunction(functionString, a + x * i) +
                  4 * calculateFunction(functionString, (2 * a + x * j) / 2) +
                  calculateFunction(functionString, a + x * (i + 1))) /
                 6);
      j += 2;
    }
  } else if (option == 2) { // simpson 3/8
    j = 1;
    y = x / 3;
    for (i = 0; i < n; i++) {
      result += (x *
                 (calculateFunction(functionString, a + x * i) +
                  3 * calculateFunction(functionString, a + y * j) +
                  3 * calculateFunction(functionString, a + y * (j + 1)) +
                  calculateFunction(functionString, a + x * (i + 1))) /
                 8);
      j += 3;
    }
  }
  return result;
}

double
trapez() { // Bu yöntemde integral n sayıda dikdörtgen kullanılarak hesaplanır.
  double result;
  int n, a, b, i;

  char functionString[128];
  printf("İntegrali hesaplanacak fonksiyonu giriniz:");
  scanf("%s", functionString);
  printf("İntegralin başlangıç değerini girin: ");
  scanf("%d", &a);
  printf("İntegralin son değerini girin: ");
  scanf("%d", &b);
  printf("n sayısını girin: ");
  scanf("%d", &n);
  double x = (double)(b - a) /
             (double)n; // int değerlerinin double olarak çevrilmesi(x=genişlik)

  if (x < 0) {
    printf("Sınırları kontrol ediniz");
  } else if (x == 0) {
    result = 0; // sınır değerler aynı ise işlem yapılmasına gerek yoktur
  } else {
    for (i = 0; i < n; i++) { // n sayıda yamuk kullanılarak
      result += (x *
                 (calculateFunction(functionString, a + x * i) +
                  calculateFunction(functionString, a + x * (i + 1))) /
                 2);
    }
  }
  return result;
}

void calculateForwardDiff(
    double y[], double diffTable[][10],
    int n) { // ileri farklar tablosunu hesaplayan fonksiyon
  int i, j;
  for (i = 0; i < n; i++) {
    diffTable[i][0] = y[i]; // ilk sütuna fonksiyonun değerleri y[]
  }
  for (j = 1; j < n; j++) {
    for (i = 0; i < n - j; i++) {
      diffTable[i][j] = diffTable[i + 1][j - 1] -
                        diffTable[i][j - 1]; // ileri farkları hesapla yerleştir
    }
  }
}

double interpolate(double x[], double diffTable[][10], int n,
                   double hesaplanacakX,
                   double h) { // x değeri için fonk değerini hesapla
  double result = diffTable[0][0];
  double u =
      (hesaplanacakX - x[0]) / h; // x değerinin h değerine göre ileri farkı
  double u_term = 1;

  for (int i = 1; i < n; i++) {
    u_term *= (u - (i - 1)) / i;
    result += u_term * diffTable[0][i]; // her terim polinomun parçası olarak
                                        // hesaplanır ve sonuca eklenir
  }

  return result;
}
float gregoryNewton() {
  int n, i;
  printf("Bilinen nokta sayısı:");
  scanf("%d", &n);
  double *xler = (double *)malloc(sizeof(double) * n);
  double *yler = (double *)malloc(sizeof(double) * n);

  for (i = 0; i < n; i++) {
    printf("%d. nokta icin x ve y degerlerini giriniz:", i + 1);
    scanf("%lf %lf", &xler[i], &yler[i]);
  }
  float x_value;
  printf("Hesaplanacak x degerini giriniz:");
  scanf("%f", &x_value);
  float h = xler[1] - xler[0];
  double diffTable[n][10];
  calculateForwardDiff(yler, diffTable, n);
  float result = interpolate(xler, diffTable, n, x_value, h);
  printf("x = %f icin y = %f\n", x_value, result);
  free(xler);
  free(yler);
  return result;
}

int main(void) {
  int choice;

  printf("Kullanmak istediğiniz yöntemi seçiniz:\n");
  printf("1.Bisection\n2.Regula Falsi\n3.Newton-Raphson\n4.Matrisin "
         "Tersi\n5.Gauss Eleminasyon\n6.Gauss Seidal\n7.Sayısal "
         "Türev\n8.Simpson\n9.Trapez\n10.Gregory Newton\n11.Çıkış\n");
  scanf("%d", &choice);
  switch (choice) {
  case 1:
    printf("Sonuc= %lf", Bisection());
    break;
  case 2:
    printf("Sonuc= %lf", RegulaFalsi());
    break;
  case 3:
    printf("Sonuc= %lf", NewtonRaphson());
    break;
  case 4:
    inverseMatrix();
    break;
  case 5:
    GaussElemination();
    break;
  case 6:
    GaussSeidal();
    break;
  case 7:
    printf("Sonuc= %lf", sayisalTurev());
    break;
  case 8:
    printf("Sonuç= %lf", Simpson());
    break;
  case 9:
    printf("Sonuç= %lf", trapez());
  case 10:
    printf("Sonuc= %lf", gregoryNewton());

  default:
    break;
  }

  return 0;
}
