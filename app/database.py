from flask import Flask, request, render_template
import os

account = {}
transaction_sig = {}


class User:
    def __init__(self, addr, arry, pk):
        self.address = addr
        self.balance = 1000
        self.mes = {}
        self.transaction_history = []
        self.need = []
        self.mulpay = []
        self.prec = arry
        self.pk = pk

    def user_page(self):
        return render_template('user.html', address=self.address, balance=self.balance, transaction_history=self.transaction_history, friend_list=self.need)

