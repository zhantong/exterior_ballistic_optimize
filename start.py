"""外弹道优化
在Windows 10(x64) Python 3.5(x64)下测试通过
"""
import math
M = [1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.50]
C = [0.383, 0.382, 0.375, 0.366, 0.356,
     0.346, 0.337, 0.328, 0.32, 0.313, 0.288]
RHO = 1.206
G = 9.8
deltat = 0.001
E = 200000


def get_tau(y):
    """气温
    由高度得到
    """
    if y <= 9300:
        return 288.9-6.328*10**-3*y
    elif y <= 12000:
        return 230.2-6.328*10**-3*(y-9300)+1.172*10**-6*(y-9300)**2
    else:
        return 221.5


def get_cs(tau):
    """声速
    由气温得到
    """
    return (1.404*287*tau)**0.5


def get_ma(v, cs):
    """马赫数
    由速度和声速得到
    """
    return v/cs


def get_cxon(ma):
    """标准弹形阻力系数
    由马赫数得到
    """
    if ma < M[0]:
        return C[0]
    if ma > M[-1]:
        return C[-1]
    for index, m in enumerate(M):
        if ma < m:
            left = ma-M[index-1]
            right = m-ma
            if left < right:
                return C[index-1]
            else:
                return C[index]


def get_i(H):
    """弹形系数
    """
    return 2.9-1.373*H+0.32*H**2-0.0267*H**3


def get_v(E, m):
    """速度
    由动能和质量得到
    """
    return (2*E/m)**0.5


def get_cx(i, cxon):
    """阻力系数
    由弹形系数和标准弹形阻力系数得到
    """
    return i*cxon


def get_S(d):
    """特征面积
    由弹径得到
    """
    return math.pi*d**2/4


def get_t(y, v, theta, i, m, S, SD):
    """迭代得到飞行时间
    """
    sum_s = 0
    count = 0
    while sum_s <= SD:
        tau = get_tau(y)
        cs = get_cs(tau)
        ma = get_ma(v, cs)
        cxon = get_cxon(ma)
        cx = get_cx(i, cxon)
        v -= (RHO*v**2/(2*m)*S*cx-G*math.sin(theta))*deltat
        theta -= G*math.cos(theta)/v*deltat
        s = v*deltat
        sum_s += s
        y += s*math.sin(theta)
        count += 1
    return deltat*count


def cal_c(m, lam_n, lam_b, lam_c):
    """迭代得到飞行时间
    即调用get_t
    """
    H = lam_n+lam_b-0.3
    i = get_i(H)
    v = get_v(E, m)
    theta = 65/360*2*3.14
    S = get_S(0.035)
    return get_t(0, v, theta, i, m, S, 2000)


def alternation():
    """获得最优解
    """
    min_t = 10000
    min_p = 10000
    m = 0.5
    for i in range(3900):
        lam_n = 0+i*0.001
        t = cal_c(0.5, lam_n, 0, 1.7)
        if t < min_t:
            min_t = t
            min_p = lam_n
    lambda_b = 0.65
    lambda_n = min_p-lambda_b
    lambda_c = 5.5-min_p
    print(m, lambda_n, lambda_b, lambda_c, min_t)
if __name__ == '__main__':
    t = cal_c(0.74, 3, 0.65, 1.7)
    print(t)
    alternation()
