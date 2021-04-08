import tkinter as tk
import tkinter.ttk as ttk

from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
from tkinter import font
import numpy as np
import pandas as pd


# Parent Classes
class MultiListbox(tk.Frame):

	def __init__(self, master, lists):
		tk.Frame.__init__(self, master)
		self.lists = []

		frame = tk.Frame(self)
		frame.pack(side=tk.RIGHT, fill=tk.Y)
		tk.Label(frame, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		sb = tk.Scrollbar(frame, orient=tk.VERTICAL, command=self._onVsb)
		sb.pack(expand=tk.YES, fill=tk.Y)

		for l, w in lists:
			frame = tk.Frame(self)
			frame.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)
			tk.Label(frame, text=l, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
			lb = tk.Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,
			                relief=tk.FLAT, exportselection=tk.FALSE, yscrollcommand=sb.set)
			lb.pack(expand=tk.YES, fill=tk.BOTH)
			self.lists.append(lb)
			lb.bind('<MouseWheel>', self._onMouseWheel)
			lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
			lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))

	def _select(self, y):
		row = self.lists[0].nearest(y)
		self.selection_clear(0, tk.END)
		self.selection_set(row)
		return 'break'

	def _button2(self, x, y):
		for l in self.lists: l.scan_mark(x, y)
		return 'break'

	def _b2motion(self, x, y):
		for l in self.lists: l.scan_dragto(x, y)
		return 'break'

	def _onVsb(self, *args):
		for l in self.lists:
			l.yview(*args)

	def _onMouseWheel(self, event):
		delta = -event.delta
		for l in self.lists:
			l.yview("scroll", delta, "units")
		return 'break'

	def curselection(self):
		return self.lists[0].curselection()

	def delete(self, first, last=None):
		for l in self.lists:
			l.delete(first, last)

	def get(self, first, last=None):
		result = []
		for l in self.lists:
			result.append(l.get(first, last))
		if last: return map([None] + result)
		return result

	def index(self, index):
		return self.lists[0].index(index)

	def insert(self, index, *elements):
		for e in elements:
			i = 0
			for l in self.lists:
				l.insert(index, e[i])
				i = i + 1

	def size(self):
		return self.lists[0].size()

	def see(self, index):
		for l in self.lists:
			l.see(index)

	def selection_anchor(self, index):
		for l in self.lists:
			l.selection_anchor(index)

	def selection_clear(self, first, last=None):
		for l in self.lists:
			l.selection_clear(first, last)

	def selection_includes(self, index):
		return self.lists[0].selection_includes(index)

	def selection_set(self, first, last=None):
		for l in self.lists:
			l.selection_set(first, last)

	def delete_all(self):
		for l in self.lists:
			l.delete(0, tk.END)

	# def sort(self, col_idx):
	# 	sorter = sorted(range(len()))


class SingleListbox(tk.Frame):

	def __init__(self, master, lists):
		tk.Frame.__init__(self, master)
		self.lists = []

		frame = tk.Frame(self)
		frame.pack(side=tk.RIGHT, fill=tk.Y)
		tk.Label(frame, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		sb = tk.Scrollbar(frame, orient=tk.VERTICAL, command=self._onVsb)
		sb.pack(expand=tk.YES, fill=tk.Y)

		frame = tk.Frame(self)
		frame.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)
		tk.Label(frame, text=lists[0], borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		lb = tk.Listbox(frame, width=lists[1], borderwidth=0, selectborderwidth=0, selectmode=tk.EXTENDED,
		                relief=tk.FLAT, exportselection=tk.FALSE, yscrollcommand=sb.set)
		lb.pack(expand=tk.YES, fill=tk.BOTH)
		self.lists.append(lb)
		lb.bind('<MouseWheel>', self._onMouseWheel)
		lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
		lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
		lb.bind('<Control-1>', lambda e, s=self: s._addselect(e.y))

	def _select(self, y):
		row = self.lists[0].nearest(y)
		self.selection_clear(0, tk.END)
		self.selection_set(row)
		return 'break'

	def _addselect(self, y):
		row = self.lists[0].nearest(y)
		self.selection_set(row)
		return 'break'

	def _button2(self, x, y):
		for l in self.lists: l.scan_mark(x, y)
		return 'break'

	def _b2motion(self, x, y):
		for l in self.lists: l.scan_dragto(x, y)
		return 'break'

	def _onVsb(self, *args):
		for l in self.lists:
			l.yview(*args)

	def _onMouseWheel(self, event):
		delta = -event.delta
		for l in self.lists:
			l.yview("scroll", delta, "units")
		return 'break'

	def curselection(self):
		return self.lists[0].curselection()

	def delete(self, first, last=None):
		for l in self.lists:
			l.delete(first, last)

	def get(self, first, last=None):
		result = []
		for l in self.lists:
			result.append(l.get(first, last))
		if last: return map([None] + result)
		return result

	def index(self, index):
		return self.lists[0].index(index)

	def insert(self, index, *elements):
		for e in elements:
			i = 0
			for l in self.lists:
				l.insert(index, e[i])
				i = i + 1

	def size(self):
		return self.lists[0].size()

	def see(self, index):
		for l in self.lists:
			l.see(index)

	def selection_anchor(self, index):
		for l in self.lists:
			l.selection_anchor(index)

	def selection_clear(self, first, last=None):
		for l in self.lists:
			l.selection_clear(first, last)

	def selection_includes(self, index):
		return self.lists[0].selection_includes(index)

	def selection_set(self, first, last=None):
		for l in self.lists:
			l.selection_set(first, last)

	def delete_all(self):
		for l in self.lists:
			l.delete(0, tk.END)

	def select_all(self):
		self.selection_set(0, tk.END)


# Derived Classes
class BinListbox(MultiListbox):
	def __init__(self, master, lists):
		tk.Frame.__init__(self, master)
		self.lists = []

		frame = tk.Frame(self)
		frame.pack(side=tk.RIGHT, fill=tk.Y)
		tk.Label(frame, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		sb = tk.Scrollbar(frame, orient=tk.VERTICAL, command=self._onVsb)
		sb.pack(expand=tk.YES, fill=tk.Y)

		for l, w in lists:
			frame = tk.Frame(self)
			frame.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)
			tk.Label(frame, text=l, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
			lb = tk.Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,
			                relief=tk.FLAT, exportselection=tk.FALSE, yscrollcommand=sb.set)
			lb.pack(expand=tk.YES, fill=tk.BOTH)
			self.lists.append(lb)
			lb.bind('<MouseWheel>', self._onMouseWheel)
			lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
			lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
			lb.bind('<Control-1>', lambda e, s=self: s._addselect(e.y))
			lb.bind('<Control-B1-Motion>', lambda e, s=self: s._addselectpressed(e.y))

	def _select(self, y):
		row = self.lists[0].nearest(y)
		self.selection_clear(0, tk.END)
		self.selection_set(row)
		self.master.highlightSelection(self.curselection())

		return 'break'

	def _addselect(self, y):
		row = self.lists[0].nearest(y)
		if row in self.curselection():
			self.selection_clear(row)
		else:
			self.selection_set(row)

		self.master.highlightSelection(self.curselection())
		return 'break'

	def _addselectpressed(self, y):
		return 'break'

	def select_all(self):
		self.selection_set(0, tk.END)

	def sort(self, indexlist):
		sorter = sorted(range(len(indexlist)), key=lambda k: indexlist[k])

		for l in self.lists:
			unsort = l.get(0, tk.END)
			sorty = [unsort[i] for i in sorter]

			l.delete(0, tk.END)
			for sortitem in sorty:
				l.insert(tk.END, sortitem)

	def itemconfig(self, item, **options):
		for l in self.lists:
			l.itemconfig(item, options)

	def itemcget(self, item, option):
		return self.lists[0].itemcget(item, option)


class GroupListbox(MultiListbox):
	def __init__(self, master, lists):
		super().__init__(master, lists)

	def _select(self, y):
		row = self.lists[0].nearest(y)
		self.selection_clear(0, tk.END)
		self.selection_set(row)
		self.master.changeBinView()
		self.master.changeMatView()
		return 'break'

	def setSize(self, size):
		idx = self.curselection()[0]
		silibo = self.lists[1]
		silibo.delete(idx)
		silibo.insert(idx, size)

		self.selection_set(idx)

	def selection_set(self, first, last=None):
		for l in self.lists:
			l.selection_set(first, last)

		self.master.changeBinView()
		self.master.changeMatView()


class MaterialListbox(MultiListbox):
	def __init__(self, master, lists):
		tk.Frame.__init__(self, master)
		self.lists = []

		frame = tk.Frame(self)
		frame.pack(side=tk.RIGHT, fill=tk.Y)
		tk.Label(frame, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
		sb = tk.Scrollbar(frame, orient=tk.VERTICAL, command=self._onVsb)
		sb.pack(expand=tk.YES, fill=tk.Y)

		for l, w in lists:
			frame = tk.Frame(self)
			frame.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)
			tk.Label(frame, text=l, borderwidth=1, relief=tk.RAISED).pack(fill=tk.X)
			lb = tk.Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,
			                relief=tk.FLAT, exportselection=tk.FALSE, yscrollcommand=sb.set)
			lb.pack(expand=tk.YES, fill=tk.BOTH)
			self.lists.append(lb)
			lb.bind('<MouseWheel>', self._onMouseWheel)
			lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
			lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))

	def _select(self, y):
		row = self.lists[0].nearest(y)
		self.selection_clear(0, tk.END)
		self.selection_set(row)
		dfidx = self.lists[0].get(row)
		self.master.changeCompoInfo(dfidx)
		return 'break'

