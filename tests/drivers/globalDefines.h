#define POINTS(number) \
  int getEvalPointsCount() {return number;} \
  extern double points[number][in_count]; \
  double getEvalPoint(int point, int col) { return points[point][col]; } \
  double points[number][in_count]

#define IN(number) \
  const int in_count = number; \
  int getInputCount() {return number;}

#define OUT(number) \
  int getOutputCount() {return number;}

int getEvalPointsCount();
double getEvalPoint(int point, int col);
int getInputCount();
int getOutputCount();

void func(NUMBER* x, NUMBER* y);
