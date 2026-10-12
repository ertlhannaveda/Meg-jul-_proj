Ötlet plusz leírás (google docsbol):

Baranyi Bartal: Forrás, Peremfeltételek & Idősoros Adatmodul
Fő feladat: A bemeneti adatok előkészítése, a kútrendszer modellezése és az időjárási adatsor integrálása.
Részfeladatok:
Helyszín kiválasztása (pl. dél-alföldi fiktív/valós kútpár: $T_{\text{in}} = 125\ ^\circ\text{C}$, m=60kg/s, $T_{\text{out}} \ge 65\ ^\circ\text{C}$).
Órás felbontású meteorológiai adatsor (külső hőmérséklet, TMY) beszerzése és feldolgozása (pl. PVGIS TMY API vagy CSV formátumban).
Hálózati betáplálási vagy helyi fogyasztási igények felmérése.
A projektterv és az adminisztráció összefogása.
Python felelősség: Az adatbetöltő és előfeldolgozó modul (data_loader.py / weather.py), amely órás idősort ad át a fő szimulációnak.

Mészáros Milán : Termodinamikai Modell & Szerves Rankine-Ciklus Számítás
Fő feladat: A fizikai alapmodell megírása névleges állapotban CoolProp segítségével.
Részfeladatok:
Az ORC (szerves Rankine-ciklus) körfolyamat állapotpontjainak leírása (elpárologtató, turbina, léghűtéses kondenzátor, tápszivattyú).
Munkaközegek előzetes vizsgálata és kiválasztása (pl. R245fa, R1233zd(E), R600a).
Komponens-hatásfokok meghatározása (izentropikus turbina- és szivattyúhatásfok, a hőcserélő minimális hőmérsékletkülönbsége).
T-S (hőmérséklet–entrópia) és p-h diagramok generálása az ellenőrzéshez.
Python: A termodinamikai alapszámítási modul (orc_thermo.py), ami adott forrás- és környezeti hőmérséklet mellett kiszámítja a teljesítményt és a fajlagos entalpiákat.


Hanna: Dinamikus Éves Szimuláció & Rendszer-optimalizálás
Fő feladat: A termodinamikai modell és az órás idősor összekapcsolása 8760 órára, valamint a működési paraméterek optimalizálása.
Részfeladatok:
A rendszer viselkedésének modellezése változó környezeti hőmérséklet mellett (nyári teljesítménycsökkenés vs. téli csúcstermelés).
Segédüzemi fogyasztások lefedése (kondenzátor ventilátorok, kút-szivattyúk).
Optimalizálás: a legkedvezőbb elpárolgási nyomás (pevap) és munkaközeg kiválasztása az éves villamosenergia-hozam maximalizálására.
Python felelősség: Az éves futtató és optimalizáló motor (annual_simulation.py / optimizer.py), éves termelési görbék és hisztogramok rajzolása.


Lelkes Barnabbás: Gazdasági, Környezeti & Társadalmi Értékelés
Fő feladat: A műszaki eredmények lefordítása pénzügyi mutatókra és fenntarthatósági mérőszámokra.
Részfeladatok:
Beruházási (CAPEX: fúrás, ORC gépészet, kiserőművi infrastruktúra) és üzemeltetési (OPEX: karbantartás, biztosítás) költségbecslés.
Gazdaságossági számítások elvégzése: Nettó Jelenérték (NPV), Belső Megtérülési Ráta (IRR), LCOE (€/MWh).
Környezeti hatáselemzés: kiváltott CO2-kibocsátás a magyar energiamixhez képest, hűtőközegek GWP vizsgálata, visszasajtolás kockázatai.
Érzékenységvizsgálat futtatása a villamosenergia-átvételi ár és a beruházási költség függvényében.
Python felelősség: A pénzügyi és környezeti kalkulátor (economics_lca.py), amely az éves energiatermelésből számol pénzáramokat és megtakarításokat.


gemini link:
https://share.gemini.google/dSxMa9TSgI5l
