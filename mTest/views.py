import time
import operator
from Bio import SeqIO
from mechanize import Browser
from bs4 import BeautifulSoup
from Bio.SeqUtils import ProtParam
from django.shortcuts import render
from django.http import HttpResponse


def homepage(request):
 return render(request, 'home.html')

def about(request):
 return render(request, 'about.html')

def tut(request):
 return render(request, 'tutorials.html')

def uselinks(request):
 return render(request, 'useful.html')

def contact(request):
 return render(request, 'contact.html')

def decide(request):
    query = request.GET['sequence']
    save(query)
    br = Browser()
    br.set_handle_robots(False)
    br.open("https://web.expasy.org/blast/")
    forms = [f.name for f in br.forms()]
    now = time.time()
    br.select_form(nr=0)
    br.form['seq'] = str(query)
    resp = br.submit()
    url = br.response()
    test = time.time()
    soup = BeautifulSoup(url, 'html.parser')
    t = soup.find('fieldset', attrs={'class': 'aln'})
    a =  str(t.find('pre'))
    f = "Identities"
    s = a.find(f)
    sub = a[s:s+27]
    inter = sub.find("(")
    intercov = sub[inter-8:inter]
    intercov = intercov.split('/')
    identity = int(intercov[0])
    AlnLength = int(intercov[1])
    cov = int((identity/AlnLength)*100)
    per = sub[inter:len(sub)]
    ide = per.replace('(', '')
    ide = ide.replace('%)', '')
    ide = int(ide)
    eta =  round((test-now)/60)
    if(cov>60 and ide>60):
        model = "Homology Modeling"
        l1 = "https://salilab.org/modeller/download_installation.html"
        l2 = "https://swissmodel.expasy.org/interactive"
        l3 = "http://www.reading.ac.uk/bioinf/IntFOLD/IntFOLD5_form.html"
        l4 = "https://www.unamur.be/sciences/biologie/urbm/bioinfo/esypred/"
        l5 = "https://modbase.compbio.ucsf.edu/modweb/"
        l6 = "http://robetta.bakerlab.org"
        link1 = "<a href={}>Modeller (Standalone)</a>".format(l1)
        link2 = "<a href={}>SwissModel</a>".format(l2)
        link3 = "<a href={}>IntFOLD</a>".format(l3)
        link4 = "<a href={}>EsyPred3D</a>".format(l4)
        link5 = "<a href={}>ModWeb</a>".format(l5)
        link6 = "<a href={}>Robetta</a>".format(l6)
        res, count, m, a, i, c, p, mc = prot()
        return render(request, 'decided.html', {'m':model, 'link1':link1, 'link2':link2, 'link3':link3, 'link4':link4, 'link5':link5, 'link6':link6, 'eta':eta, 'sequence':res, 'count': count, 'molwt':m, 'arom':a, 'instabdex':i, 'inststr':c, 'pi':p, 'mec': mc})
    elif(60>cov and cov>30 and 60>ide and ide>30):
        model = "Threading/Fold recognition Modeling"
        link1 = "a"
        link2 = "b"
        link3 = "c"
        link4 = "d"
        link5 = "e"
        link6 = "f"
        res, count, m, a, i, c, p, mc = prot()
        return render(request, 'decided.html', {'m':model,'link1':link1, 'link2':link2, 'link3':link3, 'link4':link4, 'link5':link5, 'link6':link6, 'eta':eta, 'sequence':res, 'count': count, 'molwt':m, 'arom':a, 'instabdex':i, 'inststr':c, 'pi':p, 'mec':mec})
    elif(30>cov and cov>0 and 30>ide and ide>0):
        model = "Ab Initio Modeling"
        link1 = "a"
        link2 = "b"
        link3 = "c"
        link4 = "d"
        link5 = "e"
        link6 = "f"
        res, count, m, a, i, c, p, mc = prot()
        return render(request, 'decided.html', {'m':model, 'link1':link1, 'link2':link2, 'link3':link3, 'link4':link4, 'link5':link5, 'link6':link6, 'eta':eta, 'sequence':res, 'count': count, 'molwt':m, 'arom':a, 'instabdex':i, 'inststr':c, 'pi':p, 'mec': mc})
    else:
        return render(request, 'sorry.html')
    #return render(request, 'test.html')

def save(text):
    with open("mTest/media/query.txt", "w") as f:
        f.write(text)
        f.close()

def prot():
    for seq_rec in SeqIO.parse("mTest/media/query.txt", "fasta"):
        res = str(seq_rec.seq)
    X = ProtParam.ProteinAnalysis(res)
    count = len(res)
    m = float("{0:.2f}".format(X.molecular_weight()))
    a = float("{0:.2f}".format(X.aromaticity()))
    i = float("{0:.2f}".format(X.instability_index()))
    if i > 40:
        c = "Instable"
    else:
        c = "Stable"
    p =  float("{0:.2f}".format(X.isoelectric_point()))
    mc = X.molar_extinction_coefficient()[1]
    return res, count, m, a, i, c, p, mc

'''
def webmod(request):
        seq = SeqIO.read("modelaid/media/query.txt", "fasta")
        br = Browser()                # Create a browser
        br.open("https://modbase.compbio.ucsf.edu/modweb/")
        br.select_form(nr=0)
        br.form['modweb_name'] = 'Default'
        br.form['modellerkey'] = 'MODELIRANJE'
        br.form['sequence'] = str(seq)
        resp = br.submit()
        url = br.response()
        soup = BeautifulSoup(url, 'html.parser')
        text1 = ""
        text2 = ""
        trymod = soup.find('div', attrs = {'id': 'container'})
        for row in trymod.findAll('div', attrs = {'id': 'fullpart'}):
            link = [a['href'] for a in soup.select('a[href]')][15]
            text1 = row.p.text
            text2 = row.find_all('p')[-2].text
    return render(request, 'directresult.html', {'link': link, 'text1': text1, 'text2': text2, 'sequence': seq})

def tasser(request):
    user = request.POST.get('usr_mail')
    pss = request.POST.get('pss')
    seq = SeqIO.read("modelaid/media/query.txt", "fasta")
    br = Browser()
    br.open("https://zhanglab.ccmb.med.umich.edu/I-TASSER/")
    br.select_form(nr=0)
    br.form.enctype = "multipart/form-data"
    br.form['REPLY-E-MAIL'] = id
    br.form['password'] = pss
    br.form['SEQUENCE'] = str(seq)
    resp = br.submit()
    url = br.response()
    soup = BeautifulSoup(url, 'html.parser')
    link = [a['href'] for a in soup.select('a[href]')][0]
    text = soup.find('p')
    return render(request, 'directresult.html', {'url': link, 'text': text})

def quark(request):
    user = request.POST.get('usr_mail')
    pss = request.POST.get('pss')
    seq = SeqIO.read("modelaid/media/query.txt", "fasta")
    br = Browser()
    br.open("https://zhanglab.ccmb.med.umich.edu/QUARK/")
    br.select_form(nr=0)
    br.form.enctype = "multipart/form-data"
    br.form['REPLY-E-MAIL'] = id
    br.form['PASSWORD'] = pss
    br.form['SEQUENCE'] = str(seq)
    resp = br.submit()
    url = br.response()
    text = "A confirmation email for how to retrieve your job will be send to " +id+ " in about 5 minutes"
    return render(request, 'directresult.html', {'url': u, 'text': t})
'''
