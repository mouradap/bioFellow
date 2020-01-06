from kivy.app import App
from kivy.uix.widget import Widget
from kivy.uix.button import Button
from kivy.uix.textinput import TextInput
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.uix.widget import Widget
from kivy.uix.popup import Popup
from kivy.properties import ObjectProperty
from functions import *
from kivy.uix.image import AsyncImage
from kivy.lang import Builder
from kivy.uix.dropdown import DropDown
from kivy.uix.gridlayout import GridLayout


class HomeScreen(Screen):
	pass

class WindowManager(ScreenManager):
	pass



class codingSequences(Screen):
	def cdsTranslation(self):
		cds = self.cdsSeq.text

		self.aminoacidSeq.text = translate(cds)


	def reverseTranslation(self):
		aminoacidSeq = self.aminoacidSeq.text

		self.cdsSeq.text = revTranslation(aminoacidSeq)

	def rnaDna(self):
		mRNA = self.rnaSeq.text

		self.cdsSeq.text = mRNAtoDNA(mRNA)

	def dnaRna(self):
		cds = self.cdsSeq.text

		self.rnaSeq.text = DNAtomRNA(cds)

	def clearAll(self):
		self.cdsSeq.text = ""
		self.aminoacidSeq.text = ""
		self.rnaSeq.text = ""

class orfFinding(Screen):

	def findOrf(self):
		orf = self.findOrfs.text
		orfList = orfFinder(orf)

		self.selectedOrf.text = orfList[0][1]



kv = Builder.load_file('biohub.kv')

class Credits(FloatLayout):
	pass

class BioHubApp(App):
	def build(self):
		return kv

	def creditsBtn(self):
		showCredits()


def showCredits():
	show = Credits()
	popupCredits = Popup(title = "Credits", content = show, size_hint = (.5, .5), size = (400, 400))
	popupCredits.open()


if __name__ == "__main__":
	BioHubApp().run()