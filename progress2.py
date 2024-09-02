import numpy as np
from scipy.optimize import linear_sum_assignment
import random
from enum import Enum
from datetime import datetime, timedelta
from collections import defaultdict
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading
import time

# Enum sınıfları (önceki koddan)
class AjanTipi(Enum):
    AMBULANS = 1
    ITFAIYE = 2
    POLIS = 3

class AjanDurumu(Enum):
    BOSTA = 1
    GOREVDE = 2

class Durum(Enum):
    BEKLIYOR = 1
    MUDAHALE_EDILIYOR = 2
    TAMAMLANDI = 3

class AciliyetDurumu(Enum):
    COK_ACIL = 1
    ACIL = 2
    NORMAL = 3
    AZ_ACIL = 4
    ACIL_DEGIL = 5

# Temel sınıflar (önceki koddan)
class Konum:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Ajan:
    def __init__(self, id, tip, konum):
        self.id = id
        self.tip = tip
        self.konum = konum
        self.durum = AjanDurumu.BOSTA
        self.mevcut_gorev = None

class Hasta:
    def __init__(self, id, konum, aciliyet_durumu):
        self.id = id
        self.konum = konum
        self.aciliyet_durumu = aciliyet_durumu
        self.durum = Durum.BEKLIYOR
        self.yardim_talebi_zamani = datetime.now()

class Gorev:
    def __init__(self, id, hasta, atanan_ajan):
        self.id = id
        self.hasta = hasta
        self.atanan_ajan = atanan_ajan
        self.durum = Durum.MUDAHALE_EDILIYOR
        self.tahmini_mudahale_suresi = 0

# Diğer yardımcı sınıflar (önceki koddan)
class ZamanYoneticisi:
    def __init__(self):
        self.simulasyon_zamani = datetime.now()

    def simulasyon_zamani_al(self):
        return self.simulasyon_zamani

    def bekle(self, dakika):
        self.simulasyon_zamani += timedelta(minutes=dakika)

class Optimizasyon:
    def __init__(self, simulasyon):
        self.simulasyon = simulasyon
    
    def mesafe_hesapla(self, konum1, konum2):
        return np.sqrt((konum1.x - konum2.x)**2 + (konum1.y - konum2.y)**2)
    
    def maliyet_matrisi_olustur(self, ajanlar, hastalar):
        maliyet_matrisi = np.zeros((len(ajanlar), len(hastalar)))
        for i, ajan in enumerate(ajanlar):
            for j, hasta in enumerate(hastalar):
                mesafe = self.mesafe_hesapla(ajan.konum, hasta.konum)
                aciliyet_puani = 5 - hasta.aciliyet_durumu.value
                bekleme_suresi = (self.simulasyon.zaman_yoneticisi.simulasyon_zamani_al() - hasta.yardim_talebi_zamani).total_seconds() / 60
                maliyet = mesafe * aciliyet_puani / (1 + bekleme_suresi/60)
                maliyet_matrisi[i, j] = maliyet
        return maliyet_matrisi
    
    def optimal_atama(self, ajanlar, hastalar):
        maliyet_matrisi = self.maliyet_matrisi_olustur(ajanlar, hastalar)
        satir_indeksleri, sutun_indeksleri = linear_sum_assignment(maliyet_matrisi)
        atamalar = []
        for satir, sutun in zip(satir_indeksleri, sutun_indeksleri):
            atamalar.append((ajanlar[satir], hastalar[sutun]))
        return atamalar

class LLMOnbellek:
    def __init__(self, kapasite=100):
        self.kapasite = kapasite
        self.onbellek = defaultdict(list)
    
    def anahtar_olustur(self, ajan, hastalar):
        return (ajan.id, tuple(sorted(h.id for h in hastalar)))
    
    def ekle(self, ajan, hastalar, karar):
        anahtar = self.anahtar_olustur(ajan, hastalar)
        if len(self.onbellek[anahtar]) >= self.kapasite:
            self.onbellek[anahtar].pop(0)
        self.onbellek[anahtar].append(karar)
    
    def al(self, ajan, hastalar):
        anahtar = self.anahtar_olustur(ajan, hastalar)
        if anahtar in self.onbellek and self.onbellek[anahtar]:
            return self.onbellek[anahtar][-1]
        return None

class LLMEntegrasyon:
    def prompt_olustur(self, ajan, hastalar, tum_ajanlar):
        return f"Ajan {ajan.id} için {len(hastalar)} hasta arasından seçim yap"

    def llm_cagri(self, prompt):
        return "LLM yanıtı"

    def llm_yanit_isle(self, llm_yanit):
        return {
            "secilen_hasta_id": 1,
            "gerekce": "En yakın ve en acil hasta",
            "tahmini_varis_suresi": 5
        }

class AcilDurumSimulasyonu:
    def __init__(self, harita_genislik, harita_yukseklik):
        self.harita_genislik = harita_genislik
        self.harita_yukseklik = harita_yukseklik
        self.ajanlar = []
        self.hastalar = []
        self.gorevler = []
        self.zaman_yoneticisi = ZamanYoneticisi()
        self.llm = LLMEntegrasyon()
        self.optimizasyon = Optimizasyon(self)
        self.llm_onbellek = LLMOnbellek()

    def ajan_ekle(self, tip):
        x = random.randint(0, self.harita_genislik)
        y = random.randint(0, self.harita_yukseklik)
        yeni_ajan = Ajan(len(self.ajanlar) + 1, tip, Konum(x, y))
        self.ajanlar.append(yeni_ajan)

    def hasta_ekle(self):
        x = random.randint(0, self.harita_genislik)
        y = random.randint(0, self.harita_yukseklik)
        aciliyet = random.choice(list(AciliyetDurumu))
        yeni_hasta = Hasta(len(self.hastalar) + 1, Konum(x, y), aciliyet)
        self.hastalar.append(yeni_hasta)

    def ajan_karar_al(self, ajan, mevcut_hastalar):
        onbellekten_karar = self.llm_onbellek.al(ajan, mevcut_hastalar)
        if onbellekten_karar:
            return onbellekten_karar
        
        prompt = self.llm.prompt_olustur(ajan, mevcut_hastalar, self.ajanlar)
        llm_yanit = self.llm.llm_cagri(prompt)
        karar = self.llm.llm_yanit_isle(llm_yanit)
        
        if karar:
            self.llm_onbellek.ekle(ajan, mevcut_hastalar, karar)
        
        return karar

    def gorev_olustur(self, ajan, hasta, tahmini_varis_suresi):
        yeni_gorev = Gorev(len(self.gorevler) + 1, hasta, ajan)
        yeni_gorev.tahmini_mudahale_suresi = tahmini_varis_suresi
        self.gorevler.append(yeni_gorev)
        ajan.durum = AjanDurumu.GOREVDE
        ajan.mevcut_gorev = yeni_gorev
        hasta.durum = Durum.MUDAHALE_EDILIYOR

    def simulasyonu_calistir(self, sure):
        baslangic = self.zaman_yoneticisi.simulasyon_zamani_al()
        while self.zaman_yoneticisi.simulasyon_zamani_al() - baslangic < sure:
            bosta_ajanlar = [a for a in self.ajanlar if a.durum == AjanDurumu.BOSTA]
            bekleyen_hastalar = [h for h in self.hastalar if h.durum == Durum.BEKLIYOR]
            
            if bosta_ajanlar and bekleyen_hastalar:
                atamalar = self.optimizasyon.optimal_atama(bosta_ajanlar, bekleyen_hastalar)
                
                for ajan, hasta in atamalar:
                    karar = self.ajan_karar_al(ajan, [hasta])
                    if karar:
                        self.gorev_olustur(ajan, hasta, karar["tahmini_varis_suresi"])
                        print(f"Ajan {ajan.id}, Hasta {hasta.id}'ye atandı. Gerekçe: {karar['gerekce']}")
            
            for gorev in self.gorevler:
                if gorev.durum == Durum.MUDAHALE_EDILIYOR:
                    gorev.tahmini_mudahale_suresi -= 1
                    if gorev.tahmini_mudahale_suresi <= 0:
                        gorev.durum = Durum.TAMAMLANDI
                        gorev.atanan_ajan.durum = AjanDurumu.BOSTA
                        gorev.atanan_ajan.mevcut_gorev = None
                        gorev.hasta.durum = Durum.TAMAMLANDI
                        print(f"Görev {gorev.id} tamamlandı. Ajan {gorev.atanan_ajan.id} boşta.")
            
            self.zaman_yoneticisi.bekle(1)

    def rapor_olustur(self):
        tamamlanan_gorev_sayisi = sum(1 for g in self.gorevler if g.durum == Durum.TAMAMLANDI)
        bekleyen_hasta_sayisi = sum(1 for h in self.hastalar if h.durum == Durum.BEKLIYOR)
        
        rapor = f"""Simülasyon Raporu:
        Toplam Ajan Sayısı: {len(self.ajanlar)}
        Toplam Hasta Sayısı: {len(self.hastalar)}
        Tamamlanan Görev Sayısı: {tamamlanan_gorev_sayisi}
        Bekleyen Hasta Sayısı: {bekleyen_hasta_sayisi}
        """
        return rapor

class KullaniciArayuzu:
    def __init__(self, simulasyon):
        self.simulasyon = simulasyon
        self.pencere = tk.Tk()
        self.pencere.title("Acil Durum Simülasyonu")
        self.pencere.geometry("800x600")

        self.harita_frame = tk.Frame(self.pencere)
        self.harita_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.kontrol_panel = ttk.Frame(self.pencere, padding="10")
        self.kontrol_panel.pack(side=tk.RIGHT, fill=tk.Y)

        self.harita_olustur()
        self.kontrol_paneli_olustur()

        self.guncelleme_thread = None
        self.guncelleme_devam = False

    def harita_olustur(self):
        self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.harita_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def kontrol_paneli_olustur(self):
        ttk.Button(self.kontrol_panel, text="Başlat", command=self.simulasyonu_baslat).pack(pady=5)
        ttk.Button(self.kontrol_panel, text="Durdur", command=self.simulasyonu_durdur).pack(pady=5)
        ttk.Button(self.kontrol_panel, text="Rapor Oluştur", command=self.rapor_olustur).pack(pady=5)

        ttk.Label(self.kontrol_panel, text="Simülasyon Hızı:").pack(pady=5)
        self.hiz_slider = ttk.Scale(self.kontrol_panel, from_=0.1, to=2, orient=tk.HORIZONTAL)
        self.hiz_slider.set(1)
        self.hiz_slider.pack(pady=5)

        self.durum_label = ttk.Label(self.kontrol_panel, text="Durum: Bekleniyor")
        self.durum_label.pack(pady=5)

    def haritayi_guncelle(self):
        self.ax.clear()
        self.ax.set_xlim(0, self.simulasyon.harita_genislik)
        self.ax.set_ylim(0, self.simulasyon.harita_yukseklik)

        for ajan in self.simulasyon.ajanlar:
            self.ax.plot(ajan.konum.x, ajan.konum.y, 'bo', markersize=10)

        for hasta in self.simulasyon.hastalar:
            renk = 'ro' if hasta.durum == Durum.BEKLIYOR else 'go'
            self.ax.plot(hasta.konum.x, hasta.konum.y, renk, markersize=8)

        self.canvas.draw()

    def simulasyonu_baslat(self):
        if not self.guncelleme_thread or not self.guncelleme_thread.is_alive():
            self.guncelleme_devam = True
            self.guncelleme_thread = threading.Thread(target=self.simulasyon_dongusu)
            self.guncelleme_thread.start()
            self.durum_label.config(text="Durum: Çalışıyor")

    def simulasyonu_durdur(self):
        self.guncelleme_devam = False
        if self.guncelleme_thread:
            self.guncelleme_thread.join()
        self.durum_label.config(text="Durum: Durduruldu")

    def simulasyon_dongusu(self):
        while self.guncelleme_devam:
            self.simulasyon.simulasyonu_calistir(timedelta(minutes=1))
            self.haritayi_guncelle()
            time.sleep(1 / self.hiz_slider.get())

    def rapor_olustur(self):
        rapor = self.simulasyon.rapor_olustur()
        rapor_penceresi = tk.Toplevel(self.pencere)
        rapor_penceresi.title("Simülasyon Raporu")
        rapor_text = tk.Text(rapor_penceresi, wrap=tk.WORD)
        rapor_text.insert(tk.END, rapor)
        rapor_text.pack(expand=True, fill=tk.BOTH)

    def calistir(self):
        self.pencere.mainloop()

# Kullanım örneği:
if __name__ == "__main__":
    simulasyon = AcilDurumSimulasyonu(100, 100)
    
    # Ajanları ekle
    for _ in range(5):
        simulasyon.ajan_ekle(AjanTipi.AMBULANS)
    for _ in range(3):
        simulasyon.ajan_ekle(AjanTipi.ITFAIYE)
    for _ in range(4):
        simulasyon.ajan_ekle(AjanTipi.POLIS)

    # Hastaları ekle
    for _ in range(20):
        simulasyon.hasta_ekle()

    arayuz = KullaniciArayuzu(simulasyon)
    arayuz.calistir()