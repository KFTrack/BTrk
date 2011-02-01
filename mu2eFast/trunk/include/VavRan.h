//#include <string>

class VavRan{
  public:
    static VavRan* Instance();

    float gen (float rkappa, float beta2, float ran);

  private:
    VavRan();
    VavRan(VavRan const&);             // copy constructor is private
    VavRan& operator=(VavRan const&);  // assignment operator is private
    static VavRan* m_pInstance;

    static float fninv[], u8[], v1[], v2[], v3[], v4[], v5[], v6[], v7[], v8[], w1[];
    static float edgec[], w2[], w3[], w4[], w5[], w6[], w8[], u1[], u2[], u3[], u4[], u5[], u6[], u7[];


    int i__1;
    float ret_val, r__1;

    float h__[9];
    int j, k, n;
    float s, t, v, x, y, p2, p3, q2, q3, x2, x3, y2, y3, s0, ac[14],
          hc[9], fl, fn, fu, pq, wk, xx, xy, yy, drk[5];
    int npt;
    float alfa[4], rlam, dsigm[5];
    int itype;

};
