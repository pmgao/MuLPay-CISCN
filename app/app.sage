import hashlib
import json
import time
from flask import Flask, request, render_template, flash, url_for, redirect, jsonify
import os
import database
import local_enc
import base64
import random
import string
from sm3 import *

from concurrent.futures import *
from hashlib import sha256, sha512
from typing import List
from gmpy2 import fac
from sage.stats.distributions.discrete_gaussian_polynomial import \
    DiscreteGaussianDistributionPolynomialSampler as d_gauss

class MusigL(object):
    N = 256
    q = 4197821
    Zq = GF(q)
    Rqx.< a > = PolynomialRing(Zq)
    Rzx.< b > = PolynomialRing(ZZ)
    Rrx.< c > = PolynomialRing(RR)
    Rq = Rqx.quotient(a ^ N + 1, 'x')
    R = Rzx.quotient(b ^ N + 1, 'x')
    R_ = Rrx.quotient(c ^ N + 1, 'x')

    def __init__(self, numbers):
        self.n = numbers
        self.m = 3

        self.omega = int(self.q.log(2)) + 1
        self.k = 8
        self.l = 4
        self.eta = 1

        self.sigma_b = 2 ^ (5 / 2) / sqrt(pi) * 2 ^ (2 / (self.N * self.k)) * self.N ^ (3 / 2) * sqrt(self.k * self.omega + 1)
        self.sigma_y = 2 ^ 9 / (pi * sqrt(pi)) * 2 ^ (2 / (self.N * self.k)) * self.q ^ (self.k / (self.l + self.k)) * self.N ^ 2 * sqrt(
            (self.k * self.omega + 1) * (2 + self.N + ((self.l + self.k) * self.N).log(10)))
        self.sigma_1 = self.sigma_b * (self.sigma_y) * sqrt(self.N * (2 * self.k * self.omega + 1) * (self.l + self.k))

        self.B = self.sigma_1 * sqrt(self.N * (self.l + self.k))
        self.Bn = sqrt(self.n) * self.B

        self.kappa = 23
        self.T = self.kappa ^ 2 * self.eta * sqrt(self.N * (self.l + self.k))
        self.alpha = (self.sigma_1 - 1) / self.T
        self.t = sqrt(self.N / ((pi - 1) * log(e, 2)))
        self.M = exp(pi / self.alpha ^ 2 + (pi * self.t) / self.alpha)
        self.l_ = 256
        self.A_ = self.Setup()

    def norm2(self, array):
        num = 0
        coe = []
        for poly in array:
            poly = (poly.lift()).change_ring(ZZ)
            coe.extend(poly.coefficients())
        for i in coe:
            num += int((i + self.q // 2) % self.q - self.q // 2) ^ 2
        return sqrt(num)

    def combine(seLf, a, b):
        if a < b: return 0
        return fac(int(a)) // (fac(int(b)) * fac(int(a - b)))

    def index(self, num: int, l1: int, l2: int, temp: List) -> List:
        if len(temp) == self.N:
            return temp
        if num < self.combine(l1 - 1, l2) * 2 ** l2:
            temp = [0] + temp
            return self.index(num, l1 - 1, l2, temp)
        else:
            num -= self.combine(l1 - 1, l2) * 2 ** l2
            if num < self.combine(l1 - 1, l2 - 1) * 2 ** (l2 - 1):
                temp = [1] + temp
                return self.index(num, l1 - 1, l2 - 1, temp)
            else:
                num -= self.combine(l1 - 1, l2 - 1) * 2 ** (l2 - 1)
                temp = [-1] + temp
                return self.index(num, l1 - 1, l2 - 1, temp)

    def Hagg(self, arg1, arg2):
        m = str(arg1) + str(arg2)
        temp = int(sha256(m.encode()).hexdigest(), 16)
        idx = temp % (self.combine(self.N, self.k) * 2 ** self.k)
        return self.R(self.index(idx, self.N, self.k, []))

    def Hsig(self, arg1, arg2, arg3):
        m = str(arg1) + str(arg2) + str(arg3)
        temp = int(sha512(m.encode()).hexdigest(), 16)
        idx = temp % (self.combine(self.N, self.k) * 2 ** self.k)
        return self.R(self.index(idx, self.N, self.k, []))

    def Hnon(self, arg1, arg2, arg3, r):
        m = str(arg1) + str(arg2) + str(arg3) + str(r)
        temp = int(sm3_hash(bytes_to_list(m.encode())), 16)
        temp = temp % 2 ** (self.l_ - 1)
        return bin(temp)[2:].rjust(self.l_, '0')

    def Setup(self):
        A = random_matrix(self.Rq, self.k, self.l)
        self.A_=block_matrix([[A,matrix.identity(self.k)]])
        return self.A_

    def small_poly(self, sign=None, seed = None):
        if sign:
            set_random_seed(seed)
        f = []
        for i in range(256):
            f.append(randint(-1, 1))
        return self.R(f)

    def Gen(self, seed):
        skl = [list(self.small_poly(True, seed)) for i in range(self.k + self.l)]
        sk0 = vector(self.Rq, [self.Rq(i) for i in skl])
        pk = self.A_ * sk0
        sk = vector(self.R, [self.R(i) for i in skl])
        return pk, sk

    def Agg(self, on):
        for i in range(len(on)):
            if not on[i][0]:
                return None
        z = sum([on[i][0] for i in range(len(on))])
        w = on[0][1]
        z_ = []
        for i in z:
            z_.append(self.Rq(i.lift()))
        return w, vector(self.Rq, z_)

    def KAgg(self, L):
        t_ = 0
        for i in range(len(L)):
            t_ += self.Rq(self.Hagg(self.order(L), L[i]).lift()) * L[i]
        return t_

    def Ver(self, pk, sigma, miu):
        w_, z = sigma
        t_ = pk
        c = self.Rq(self.Hsig(w_, miu, t_).lift())
        if self.A_ * z - c * t_ == w_ and self.norm2(list(z)) <= self.Bn:
            return True
        else:
            return False

    def gauss_v(self, sigma, length):
        def arr_append(j):
            v[j] = d()
            v_[j] = list(v[j])

        d = d_gauss(self.R, self.N, sigma)
        v = [None] * length
        v_ = [None] * length

        pool = ThreadPoolExecutor(max_workers=4)
        all_task = []
        for i in range(length):
            all_task.append(pool.submit(arr_append, i))
        wait(all_task, return_when=ALL_COMPLETED)
        pool.shutdown()

        return vector(self.R, v), v_

    def Samp(self, r, choice):
        set_random_seed(r)
        if choice:
            p = d_gauss(self.Rq, self.N, self.sigma_b)
        else:
            p = d_gauss(self.R, self.N, self.sigma_b)
        return p()

    def RejSamp(self, v, z, b):
        def ComputeN1():
            N1=0
            for row in sz:
                f = row[0].lift()
                arr = f.coefficients()
                brr = []
                for i in range(len(arr)):
                    brr.append((arr[i] % self.q) ** 2)
                N1 += sum(brr)
            return N1
        def ComputeN2():
            N2=0
            for row in szc:
                f = row[0].lift()
                arr = f.coefficients()
                brr = []
                for i in range(len(arr)):
                    brr.append((arr[i]) ** 2)
                N2 += sum(brr).lift()
            return N2

        Sigma = 0
        for ele in b:
            Sigma += ele ^ 2

        Sigma = round(self.sigma_1 ^ 2) + Sigma * round(self.sigma_y ^ 2)
        s_Sigma = matrix.identity(self.l + self.k)
        s_Sigma = matrix(self.Rq, s_Sigma)
        s_Sigma[0, 0] = Sigma ^ (-1)

        h_Sigma = round(self.sigma_1 ^ 2) * matrix.identity(self.l + self.k)
        sh_Sigma = h_Sigma.cholesky()

        z_ = matrix(self.Rrx, z).T
        z_ = matrix(self.R_, z_)
        sz = matrix(self.R_, sh_Sigma.inverse() * z_)

        tpl = list(z - v)
        tpv = []
        for i in tpl:
            tpv.append(self.Rq(i.lift()))
        temp = matrix(self.Rq, tpv).T
        szc = matrix(self.Rq, s_Sigma * temp)

        pool = ThreadPoolExecutor(max_workers=2)
        all_task=[pool.submit(ComputeN1),pool.submit(ComputeN2)]

        set_random_seed()
        X = RealDistribution('uniform', [0, 1])
        lnrho = ln(X.get_random_element())

        wait(all_task, return_when=ALL_COMPLETED)
        N1=all_task[0].result()
        N2=all_task[1].result()
        pool.shutdown()

        lnPr = -pi * (N1 - N2) - ln(self.M)
        if lnrho >= min(0, lnPr):
            return 0
        return 1

    def SignOff(self, A_, t):
        def arr_append(j):
            temp_arr[j] = self.gauss_v(self.sigma_y, self.l + self.k)
            y[j] = temp_arr[j][0]
            y_[j] = vector(self.Rq, [self.Rq(item) for item in temp_arr[j][1]])

        temp_arr = [None] * (self.m)
        y = [None] * (self.m)
        y_ = [None] * (self.m)

        pool = ThreadPoolExecutor(max_workers=self.m)
        all_task = []
        for i in range(1, self.m + 1):
            all_task.append(pool.submit(arr_append, i))

        tp1, tp2 = self.gauss_v(self.sigma_1, self.l + self.k)
        y[0] = tp1
        y_[0] = vector(self.Rq, [self.Rq(item) for item in tp2])
        wait(all_task, return_when=ALL_COMPLETED)
        pool.shutdown()

        com = []
        for i in range(len(y)):
            com.append(A_ * y_[i])
        off = [t, com]
        st = []
        st.extend(y)
        st.append(com)
        return off, st

    def order(self, mat):
        mat = matrix(mat)
        n_rows = mat.nrows()
        n_cols = mat.ncols()

        mat_copy = [[0 for col in range(n_cols)] for row in range(n_rows)]

        arr = []
        for i in range(n_rows):
            weight = 0
            for j in range(n_cols):
                f = mat[i][j].lift().change_ring(ZZ)
                weight += self.q ^ (j * self.N) * f(self.q)
            arr.append((weight, i))
        for i in range(n_rows):
            arr.sort(key=lambda x: x[0])
        for i in range(n_rows):
            row = arr[i][1]
            mat_copy[i] = mat[row]

        return mat_copy

    # msgs -> off
    def SignOn(self, idx, st, msgs, sk, miu, all_pk, rd):
        def assist1():
            for k in range(self.n):
                W.append(t_array[k])
                W.extend(com_array[k])

        def assist2():
            def assist_sample1():
                b1.extend([self.Samp(hash_res, True)] * (self.m - 1))
            def assist_sample2():
                b2.extend([self.Samp(hash_res, False)] * (self.m - 1))

            hash_res = self.Hnon(self.order(list(W)), miu, t_, rd)
            pool_sample = ThreadPoolExecutor(max_workers=4)
            all_tasks = [pool_sample.submit(assist_sample1), pool_sample.submit(assist_sample2)]
            wait(all_tasks, return_when=ALL_COMPLETED)
            pool_sample.shutdown()

        def assist3():
            for k in range(self.m):
                temp = 0
                for j in range(self.n):
                    temp += com_array[j][k]
                w.append(temp)

        t_array = []
        com_array = []
        for i in range(len(msgs)):
            t_array.append(msgs[i][0])
            com_array.append(msgs[i][1])
        for i in range(len(t_array)):
            if t_array[i] != all_pk[i] or (t_array[i] == t_array[idx] and i != idx):
                print("[-] fail, pk is error!")
                return None
        t_ = self.KAgg(t_array)

        w = []
        b1 = [self.Rq(1)]
        b2 = [self.R(1)]
        W = []

        ppool = ThreadPoolExecutor(max_workers=4)
        all_task = [ppool.submit(assist1),ppool.submit(assist3)]
        a = self.R(self.Hagg(self.order(t_array), t_array[idx]).lift())
        wait(all_task, return_when=ALL_COMPLETED)
        ppool.shutdown()
        assist2()

        for item in w[-1]:
            try:
                item ^ (-1)
            except:
                print("[-] fail, no inverse")
                return None
        w_ = 0
        y_ = 0
        for i in range(self.m):
            w_ += b1[i] * w[i]
            y_ += b2[i] * st[i]
        c = self.R(self.Hsig(w_, miu, t_).lift())
        v = c * a * sk
        z = v + y_

        #if self.RejSamp(v, z, b1) == 0:
        #    print("[-] fail, invalid sample")
        #    return None
        on = (z, w_)

        return on

    def gen_key(self):
        self.A_ = self.Setup()
        pk, sk = [], []
        # Generate pk, sk for all object
        for i in range(self.n):
            pki, ski = self.Gen(self.A_)
            pk.append(pki)
            sk.append(ski)
        return self.A_, pk, sk


    def Offline_assist(self,A_, j, sk, pk):
        offi, sti = self.SignOff(A_, sk[j], pk[j])
        print("[+] Off Siganature {} generate successfully".format(j))
        return offi, sti

    def Offline(self, A_, sk, pk):
        off, st = [], []
        pools = ProcessPoolExecutor(max_workers=self.n)
        all_task = []
        for i in range(self.n):
            all_task.append(pools.submit(self.Offline_assist, A_, i, sk, pk))
        wait(all_task, return_when=ALL_COMPLETED)

        for ele in all_task:
            off.append(ele.result()[0])
            st.append(ele.result()[1])
        pools.shutdown()
        return off, st

    def Online_assist(self,j,st,off,sk,m,pk):
        temp = self.SignOn(j, st, off, sk, m, pk)
        if not temp:
            print("[-] Siganature {} generate fail".format(j))
            return None
        print("[+] Siganature {} generate successfully".format(j))
        return temp

    def Online(self, st, off, sk, m, pk):
        on = []
        pools = ProcessPoolExecutor(max_workers=self.n)
        all_task = []
        for i in range(self.n):
            all_task.append(pools.submit(self.Online_assist, i, st, off, sk, m, pk))
        wait(all_task, return_when=ALL_COMPLETED)

        for ele in all_task:
            on.append(ele.result())
        pools.shutdown()
        return on

    def verify(self, on, pk, msg):
        # Calulate aggregate signature
        sigma = self.Agg(on)
        # Calculate aggregate public key
        agg_pk = self.KAgg(pk)
        # Verify
        return self.Ver(agg_pk, sigma, msg)

    def signature(self):
        while True:
            off, st = self.Offline(self.A_, self.sk, self.pk)
            on = self.Online(st, off, self.sk, self.message, self.pk)
            if on:
                print("[+] Generate all signature")
                break
        print("Verify signature result: {}".format(self.verify(on, self.pk)))

def auth(pk, sk, system):
    ### finish a challenge ###
    c = system.Rq(sys.Hsig(1, 2, 3).lift())
    r = d_gauss(system.Rq, system.N, system.sigma_y)
    y = []
    for _ in range(system.k+system.l):
        y.append(r())
    y = vector(system.Rq, y)
    pk_ = []
    sk_ = []
    for item in pk:
        pk_.append(system.Rq(item))
    for item in sk:
        sk_.append(system.Rq(item))
    pk_ = vector(system.Rq, pk_)
    sk_ = vector(system.Rq, sk_)
    sk_ = c*sk_+y
    t = system.A_*y
    print("debug2", system.A_*y-t)
    return system.A_*sk_ == c*pk_+t


def k2list(f):
    tp = []
    for item in f:
        tp.append(list(item))
    return tp


table = string.ascii_uppercase+string.ascii_lowercase
def random_word():
    return ''.join(random.sample(table, k=5))

check_t = '0123456789abcdef'
def check(id):
    tp1 = id[:2]
    id = id[2:]
    if tp1 != '0x':
        return False
    if all([id[i] in check_t for i in range(len(id))]):
        return True
    else:
        return False

def check_sign(item, trc, h):
    all_pk = trc['all_pk']
    msg = item['msg']
    if item['time'] == len(item['sign']):
        if sys.verify(item['sign'], all_pk, msg):
            acc = database.account[trc['source_addr']]
            acc.balance = round(acc.balance - trc['amount'], 5)
            maxl = max(len(usr.address)+29,len(str(trc['amount']))+len(trc['target_addr'])+38)

            out = f"Successful transaction\n"
            out += f"{'#' * (maxl+1)}\n"
            out += f"# transaction: {usr.address} sent {str(trc['amount'])}".ljust(maxl)+f"#\n"
            out += f"# ETH to {trc['target_addr']} has on the chain ".ljust(maxl)+f"#\n"
            out += f"{'#' * (maxl+1)}"

            print(out)
            acc.transaction_history.append({"type": "Sent", "amount": "{} ETH to {}".format(trc['amount'], trc['target_addr'])})
            database.transaction_sig.pop(h, None)
            for _ in acc.need:
                try:
                    database.account[_].mulpay.remove(trc)
                except:
                    continue




sys = MusigL(2)
usr = None
pk = None
sk = None
usr_st = None
word = None
app = Flask(__name__, template_folder='./static/page')
app.config["SECRET_KEY"] = os.urandom(32).hex()


@app.route("/")
def index():
    try:
        f = open('account.bin', 'rb')
        f.close()
        return redirect(url_for('login'))
    except:
        return redirect(url_for('welcome'))


@app.route("/welcome")
def welcome():
    return render_template('new.html')


@app.route("/register", methods=['POST', 'GET'])
def register():
    if request.method == 'GET':
        return render_template('NewWallet.html')
    global word
    word = [random_word() for _ in range(12)]
    pwd = request.form.get("password")
    pk, sk = sys.Gen(int(sha256(' '.join(word).encode()).hexdigest(), 16))
    off, st = sys.SignOff(sys.A_, pk)
    st1 = []
    for item in st[:-1]:
        st1.append(k2list(item))
    st2 = []
    for item in st[-1]:
        st2.append(k2list(item))
    enc_st = local_enc.enc_key(str(st1)+"+"+str(st2), str(k2list(sk)))
    address = '0x'+hashlib.shake_128(str(k2list(pk)).encode()).hexdigest(int(10))
    obj = database.User(address, (off, enc_st), pk)
    database.account[address] = obj
    f = open("account.bin", mode='a')
    f.close()
    f = open("account.bin", 'wb')
    f.write(base64.b64encode(local_enc.enc_key(str(k2list(pk))+"+"+str(k2list(sk)), pwd)))
    f.close()
    return redirect(url_for('show'))


@app.route("/show")
def show():
    return render_template('NewWallet2.html', wd=word)


@app.route("/login", methods=['POST', 'GET'])
def login():
    if request.method == 'GET':
        return render_template('login.html')
    password = request.form.get("pwd")
    enc_account = base64.b64decode(open('account.bin', 'rb').read())
    try:
        account = local_enc.dec_key(enc_account, password).decode().split('+')
        pk_, sk_ = account
        address = '0x'+hashlib.shake_128(pk_.encode()).hexdigest(int(10))
        global pk, sk
        pk = eval(pk_)
        sk = eval(sk_)
        if auth(pk, sk, sys):
            global usr
            usr = database.account[address]
            print("#################################################")
            print(f'# user {usr.address} login successful  #')
            print("#################################################")
            st1, st2 = local_enc.dec_key(usr.prec[1], str(sk)).decode().split('+')
            st1 = eval(st1)
            st2 = eval(st2)
            global usr_st
            usr_st = []
            for item in st1:
                st1_ = []
                for pl in item:
                    st1_.append(sys.R(pl))
                usr_st.append(vector(sys.R, st1_))
            tp = []
            for item in st2:
                st2_ = []
                for pl in item:
                    st2_.append(sys.Rq(pl))
                tp.append(vector(sys.Rq, st2_))
            usr_st.append(tp)
            return redirect(url_for('usr_state'))
        else:
            return render_template('login.html', mes="Incorrect password")
    except Exception as e:
        print(str(e))
        return render_template('login.html', mes="Incorrect password")


@app.route("/pay", methods=['POST', 'GET'])
def pay_state():
    if request.method == 'GET':
        return render_template('prev_pay.html')
    data = request.get_data()
    data_dict = {}
    for item in data.decode().split('&'):
        key, value = item.split('=')
        data_dict[key] = value
    total_amount = float(data_dict['totalAmount'])
    address = data_dict['address']
    global usr
    if usr.balance >= total_amount:
        if len(usr.need):
            off = [database.account[item].prec[0] for item in usr.need]
            st = usr_st
            all_pk = [database.account[item].pk for item in usr.need]
            trc = {"source_addr": usr.address, "target_addr": address, "amount": total_amount,\
                "token": sha256(str(time.time()).encode()).hexdigest(), "off": off, "all_pk": all_pk}
            h = sm3_hash(bytes_to_list(str(trc).encode()))
            database.transaction_sig[h] = {"time": len(usr.need), "sign": [], "msg": "aaaa"}
            for item in usr.need:
                database.account[item].mulpay.append(trc)
            return redirect('usr_state')
        else:
            usr.balance = round(usr.balance - total_amount, 5)
            usr.transaction_history.append({"type": "Sent", "amount": "{} ETH to {}".format(total_amount, address)})
            maxl = max(len(usr.address)+29,len(str(total_amount))+len(address)+38)

            out = f"Successful transaction\n"
            out += f"{'#' * (maxl+1)}\n"
            out += f"# transaction: {usr.address} sent {str(total_amount)}".ljust(maxl)+f"#\n"
            out += f"# ETH to {address} has on the chain ".ljust(maxl)+f"#\n"
            out += f"{'#' * (maxl+1)}"

            print(out)
            return redirect('usr_state')
    else:
        return redirect('login')


@app.route("/usr_state", methods=['POST', 'GET'])
def usr_state():
    global usr
    if request.method == 'GET':
        maxl = max(len(usr.address), len(str(usr.balance)))
        us = f'user state:\n'
        us += f'{"#"*(maxl+14)}\n'
        us += f'# address: {usr.address.ljust(maxl)}  #\n'
        us += f'# balance: ${str(usr.balance).ljust(maxl)} #\n'
        us += f'{"#"*(maxl+14)}\n'
        print(us)
        print("Waiting transaction info: ")
        print(database.transaction_sig)
        return render_template('user.html', address=usr.address, balance=usr.balance,\
                               transaction_history=usr.transaction_history, friend_list=usr.need, p_list = usr.mulpay)
    return render_template('user.html', address=usr.address, balance=usr.balance,\
                           transaction_history=usr.transaction_history, friend_list=usr.need, p_list = usr.mulpay)


@app.route("/forget", methods=['POST', 'GET'])
def forget():
    if request.method == 'GET':
        return render_template('forget.html')
    password = request.form.get("pwd1")
    com_password = request.form.get("pwd2")
    wd = request.form.get("word")
    try:
        assert password == com_password
        global word
        word = wd.split(' ')
        pk, sk = sys.Gen(int(sha256(' '.join(word).encode()).hexdigest(), 16))
        f = open("account.bin", mode='a')
        f.close()
        f = open("account.bin", 'wb')
        f.write(base64.b64encode(local_enc.enc_key(str(k2list(pk)) + "+" + str(k2list(sk)), password)))
        f.close()
        return redirect(url_for('login'))
    except:
        return render_template('forget.html', mes="Password confirmation does not match")


@app.route('/add_friend', methods=['POST', 'GET'])
def add_friend():
    friend_id = request.form['friend_id']
    if friend_id not in usr.need and check(friend_id):
        usr.need.append(friend_id)
    return redirect('usr_state')


@app.route('/delete_friend', methods=['POST', 'GET'])
def delete_friend():
    friend_id = request.form['friend_id']
    if friend_id in usr.need and check(friend_id):
        usr.need.remove(friend_id)
    return redirect('usr_state')


@app.route('/mulsign', methods=['POST', 'GET'])
def mulsign():
    token = request.form['token']
    for item in usr.mulpay:
        if item['token'] == token:
            global sk
            sk_ = []
            for pl in sk:
                sk_.append(sys.R(pl))
            sk_ = vector(sys.R, sk_)
            off = item['off']
            st = usr_st
            all_pk = item['all_pk']
            h = sm3_hash(bytes_to_list(str(item).encode()))
            msg = database.transaction_sig[h]['msg']
            on = sys.SignOn(off.index(usr.prec[0]), st, off, sk_, msg, all_pk, 0)
            payment = database.transaction_sig[h]
            payment['sign'].append(on)
            check_sign(payment, item, h)
    return redirect('usr_state')

app.run()

