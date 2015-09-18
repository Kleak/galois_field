library galois_field;

import "dart:typed_data";

List<int> GF_EXP = new List.filled(512, 1);
List<int> GF_LOG = new List.filled(256, 0);

int gfDivide(int x, int y) {
  if (y == 0) {
    throw "Divide by 0";
  }
  if (x == 0) {
    return 0;
  }
  return GF_EXP[GF_LOG[x] + 255 - GF_LOG[y]];
}

int gfInverse(int y) {
  return gfDivide(1, y);
}

int gfMultiply(int x, int y) {
  if (x == 0 || y == 0) return 0;
  return GF_EXP[GF_LOG[x] + GF_LOG[y]];
}

/**
 * Addition two polynomials (using exclusive-or, as usual)
 */
List<int> gfPolynomialAdd(List<int> p, List<int> q) {
  List<int> r = new List.filled(p.length > q.length ? p.length : q.length, 0);
  for (int i = 0; i < p.length; i++) {
    r[i + r.length - p.length] = p[i];
  }
  for (int i = 0; i < q.length; i++) {
    r[i + r.length - q.length] ^= q[i];
  }
  return <int>[]..addAll(r);
}

/**
 * Fast polynomial division by using Extended Synthetic Division and optimized for GF(2^p) computations
 * (doesn't work with standard polynomials outside of this galois field, see the Wikipedia article for generic algorithm).
 */
List<int> gfPolynomialDivide(List<int> dividend, List<int> divisor) {
  Uint8List msg_out = new Uint8List.fromList(dividend);
  for (int i = 0; i < dividend.length - divisor.length - 1; i++) {
    int coef = msg_out[i];
    if (coef != 0) {
      for (int j = 1; j < divisor.length; j++) {
        msg_out[i + j] += -divisor[j] * coef;
      }
    }
  }
  int separator = divisor.length - 1;
  return msg_out.sublist(msg_out.length - separator);
}

/**
 * Evaluate a polynomial at a particular value of x, producing a scalar result
 */
int gfPolynomialEval(List<int> p, int x) {
  int y = p[0];
  for (int i = 1; i < p.length; i++) {
    y = gfMultiply(y, x) ^ p[i];
  }
  return y;
}

/**
 * Multiplies two polynomials
 */
List<int> gfPolynomialMultiply(List<int> p, List<int> q) {
  List<int> r = new List.filled(p.length + q.length - 1, 0);
  for (int j = 0; j < q.length; j++) {
    for (int i = 0; i < p.length; i++) {
      r[i + j] ^= gfMultiply(p[i], q[j]);
    }
  }
  return <int>[]..addAll(r);
}

/**
 * Multiplies a polynomial by a scalar
 */
List<int> gfPolynomialScale(List<int> p, int x) {
  List<int> r = new List.filled(p.length, 0);
  for (int i = 0; i < p.length; i++) {
    r[i] = gfMultiply(p[i], x);
  }
  return <int>[]..addAll(r);
}

void initTables([int size = 255]) => _initTables(0x11d, size);

/**
 * Precompute the logarithm and anti-log tables for faster computation later, using the provided primitive polynomial
 */
void _initTables(int prim, int size) {
  int x = 1;
  for (int i = 1; i < size; i++) {
    x <<= 1;
    if (x & (size + 1) == size + 1) {
      x ^= prim;
    }
    GF_EXP[i] = x;
    GF_LOG[x] = i;
  }
  for (int i = size; i < (size + 1) * 2; i++) {
    GF_EXP[i] = GF_EXP[i - size];
  }
}
