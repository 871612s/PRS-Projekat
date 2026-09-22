## Performanse računarskih sistema

## Domaći zadatak za školsku 2025/2026. godinu

Predmetni profesor: dr Jelica Protić

Predmetni asistenti: Predrag Obradović i Matija Dodović

## Uvod

analize performansi jednog jednostavnog sistema. Preporučujemo da se njihovom rešavanju pristupi redom kojim su navedeni. Smatramo da će na taj način implementirani sistem biti logično strukturiran, a proces izrade projekta fokusiran.

Prva kategorija bavi se analitičkim rešavanjem problema koristeći teorijska razmatranja

obrađena na predmetu. Kao i uvek u stohastičkim proračunima, sistem se idealizuje i rezultati koji se dobijaju su relativno dobra ocena očekivanih vrednosti, ali bez merenja ili simulacije ne možemo biti sigurni koliko dobro opisuju sistem koji razmatramo jer ne znamo da li važe stohastičke pretpostavke.

Zahtevi iz druge kategorije vezani su za simulaciju sistema. Umesto idealizovanih

stohastičkih pretpostavki, svaki servisni centar se nezavisno modeluje. Simulacioni metod je primenljiv i u rešavanju sistema koji nisu rešivi analitički.

Treća kategorija sadrži stavke vezane za dokumentovanje, vizuelizaciju, upoređivanje i

tumačenje rezultata dobijenih analitičkim rešavanjem i simulacijom. Istinsko razumevanje rada realnog sistema može proisteći samo iz detaljnog upoređivanja analitičke i simulacione metode, pa se i u ovom veoma jednostavnom problemu tome pridaje veliki značaj.

Zahtevi u okviru ovog projekta grupisani su u 3 kategorije koje odražavaju različita lica

## Opšte napomene

Program je dozvoljeno pisati u jednom od sledećih programskih jezika: Python, C, C++, Java,

C#. Domaći zadatak se radi i brani samostalno. Ukoliko se na odbrani utvrdi da student nije samostalno radio zadatak, dobiće -5 poena koji se dodaju na konačni broj poena dobijenih na osnovu definisane formule za ocenjivanje i važe dva ispitna roka nakon odbrane (uključujući ispitni rok u kojem je bila odbrana domaćeg).

Arhivu sa izvornim kodom programa i prevedenim programom, spremnim za izvršavanje,

rezultatima simulacije i analize i traženim dijagramima, kao i dokumentovanom analizom rešenja potrebno je poslati preko odgovarajuće veb forme, koja će biti blagovremeno aktivirana. Štampanu dokumentaciju je potrebno predati prilikom odbrane.

Rok za predaju domaćeg zadatka objavljuje se putem mejling liste. Posle toga, forma za slanje

domaćih se gasi. Datum odbrane će biti određen naknadno, a u zavisnosti od rasporeda ispita u ispitnom roku.


## Opis sistema

Multiprogramski računar modeluje se otvorenom mrežom koja radi u stacionarnom režimu.

Sadrži procesor, tri sistemska diska i 𝐾 korisničkih diskova. Poasonov tok poslova intenziteta 𝛼 pristiže na procesor. Procesorska obrada traje u proseku 𝑆𝑝 = 4𝑚𝑠, prosečno vreme opsluživanja za sistemske diskove iznosi po 𝑆𝑑1 = 10𝑚𝑠, 𝑆𝑑2 = 𝑆𝑑3 = 15𝑚𝑠, a za korisničke diskove po 𝑆𝑑𝑘 = 25𝑚𝑠. Sva vremena imaju eksponencijalnu raspodelu.

Posle procesorske obrade u 20% slučajeva zahteva se pristup prvom sistemskom disku, u 15%

slučajeva drugom, u 10% slučajeva trećem, a u 40% slučajeva pristupa se nekom od 𝐾 korisničnih diskova (sa jednakom verovatnoćom za svaki disk). U 15% slučajeva posle obrade na procesoru posao se vraća u red za čekanje procesora.

Posle pristupa nekom sistemskom disku, u 35% slučajeva ne vrši se nikakva obrada već se

procesi ponovo vraćaju u procesorski red. U 25% slučajeva ponavlja se obrada na tom sistemskom disku. U preostalih 40% slučajeva pristupa se nekom od 𝐾 korisničkih diskova (sa jednakom verovatnoćom za svaki disk).

Posle pristupa korisničkom disku proces napušta sistem.

## Zadaci

resursa, protoke kroz resurse, prosečan broj poslova u svakom od resursa ovog sistema i vreme odziva ovog sistema sa centralnim serverom za K od 2 do 5. Ovi parametri se određuju za različite vrednosti ulaznog toka 𝑟 ⋅ 𝛼𝑚𝑎𝑥, 𝑟 ∈ {0.30, 0.55, 0.80, 1.00}. Odrediti kritični resurs u sistemu.

Za otvorenu mrežu kojom se modelira ovaj računarski sistem, potrebno je odrediti iskorišćenja

## Analitičko rešavanje (10p)

- 1) (4p) Obezbediti funkcionalnost koja na osnovu definisanih ulaznih parametara sistema (topologije sistema date kroz matricu verovatnoća tranzicija P i vektora brzina servera) matričnom metodom za rešavanje otvorenih mreža prikazanom na predavanjima određuje protoke kroz servere u funkciji od 𝛼. Izračunavanje protoka sprovesti za svako K (broj korisničkih diskova) i odnose protoka i intenziteta ulaznog toka upisati u prvi izlazni fajl (nazvati fajl logično, recimo “protoci_analiticki”).

- 2) (3p) Odrediti granične vrednosti intenziteta ulaznog toka 𝛼𝑚𝑎𝑥 za koje je sistem u stacionarnom režimu (u zavisnosti od K). Za svaku vrednost K iz datog opsega odrediti kritični resurs i dati grafik zavisnosti 𝛼𝑚𝑎𝑥(𝐾).

- 3) (3p) Korišćenjem Džeksonove teoreme odrediti sve tražene parametre i upisati ih u drugi izlazni fajl (nazvati fajl logično, recimo “rezultati_analiticki”). Analizu sprovesti za svako K (broj korisničkih diskova).


## Simulacija (10p)

Potrebno je implementirati funkcionalnost koja će simulirati dati sistem i na osnovu rezultata

simulacije odrediti i prikazati (upisati u treći izlazni fajl, nazvan logično, recimo “rezultati_simulacija”) sve tražene parametre sistema za svako K (broj korisničkih diskova). Osim variranja K, varira se i intenzitet ulaznog toka. Intenziteti za koje treba izvršiti simuliranje su dati sa 𝑟 ⋅ 𝛼𝑚𝑎𝑥, 𝑟 ∈ {0.30, 0.55, 0.80, 1.00}.

Potrebno je obezbediti mogućnost podešavanja simuliranog vremena rada sistema (u

minutima), a podrazumevano vreme trajanja simulacije je 0.5h (ne realnog vremena rada programa, već simuliranog vremena rada sistema!).

Simulaciju treba ponoviti 100 puta za svaku od vrednosti parametara r i K i rezultate

usrednjiti, pa usrednjene rezultate upisati u četvrti fajl (nazvan logično, recimo “rezultati_simulacija_usrednjeno”). Naravno, usrednjavaju se rezultati simulacije za isto r i K!

## Analiza i dokumentovanje rezultata (10p)

Potrebno je napisati dokumentaciju koja detaljno objašnjava metod simulacije i analitičko

rešavanje problema i upoređuje rezultate dobijene analitičkom metodom, simulacijom i usrednjene rezultate više simulacija.

U dokumentaciji priložiti tabelarne izveštaje relativnog odstupanja traženih parametara

sistema dobijenih analičkom metodom od rezultata simulacije i od usrednjenih rezultata više simulacija. Tabelarni izveštaj treba da prikazuje vrednosti i relativna odstupanja traženih parametara u zavisnosti od broja korisničkih diskova u sistemu za vrednosti parametra r date u opisu simulacije. Šta se može zaključiti o rezultatima simulacije i usrednjenim rezultatima više simulacija? Koji od njih imaju manje relativno odstupanje od rezultata dobijenih analitičkom metodom? Kako i zašto broj izvršenih simulacija utiče na relativno odstupanje?

Dokumentacija treba da sadrži i sledeće dijagrame, konstruisane na osnovu rezultata analitičkog rešavanja:

- a) Dijagrame zavisnosti iskorišćenja procesora, korisničkog diska i sistemskog diska od K (crtati različitim bojama i simbolima na istom grafiku i priložiti legendu koja jasno govori šta koja kriva predstavlja) za vrednosti parametra r date u opisu simulacije (ukupno 4 grafika, po jedan za svaku vrednost r). Kao i do sada, K se kreće od 2 do 5.

- b) Dijagrame zavisnosti vremena odziva procesora, korisničkog diska i sistemskog diska od K (crtati različitim bojama i simbolima na istom grafiku i priložiti legendu koja jasno govori šta koja kriva predstavlja) za vrednosti parametra r date u opisu simulacije (ukupno 4 grafika, po jedan za svaku vrednost r). Kao i do sada, K se kreće od 2 do 5.

- c) Dijagram zavisnosti vremena odziva sistema od K za vrednosti parametra r date u opisu simulacije (ukupno 4 grafika, po jedan za svaku vrednost r). Kao i do sada, K se kreće od 2 do 5.

Odgovoriti i na pitanje koji je resurs kritičan za svaku kombinaciju parametara r i K.


Datum izmene

## Zapisnik izmena

Izmena
