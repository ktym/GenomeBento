# 未来のゲノム弁当食材の追加のしかた

## ToyToyさんを救え！

「この食材のゲノムはまだ決まっていないのかな？」と思ったら、

[ゲノム弁当 - Google スプレッドシート](https://docs.google.com/spreadsheets/d/1iyh41eZSsOWCkZOS0b5iqjZM4jA6osboCqXFHqVtLx8/edit#gid=335143120)

を開いて探してみてください。

載ってなかったからといって、本当にまだ世界中で誰もその食材のゲノムを読んでいないのかは必ずしもわかりません。

### Wikipediaなどで食材の学名を検索

今回は、例として（すでにゲノムが決まっていますが）山口の食材として有名な「トラフグ」を検索してみます。
右上のインフォボックスの中から生物種の学名を見つけます。

![トラフグ - Wikipedia](images/Wikipedia-Trafugu.png)

トラフグの学名は Takifugu rubripes であることが分かりました。

### 生物系統分類データベースを検索

[NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) データベースで生物種 ID と系統分類を学名から調べることができます。

![Taxonomy Browser](images/Taxonomy-Torafugu.png)

これでトラフグのタクソノミー（系統分類）の ID が 31033 であることが分かります。
また、ゲノムへのリンクから下記の NCBI Genome データベースに移動して、ゲノム配列を取得することができます。

### ゲノムデータベースを検索

ゲノムのデータベースは、いろいろな生物種がまとめられたもの、ある生物種について詳しくまとめられたものなど、世界中に数多くあります。

* 解読されたゲノム配列を網羅的に収集したデータベース [NCBI Genome](https://www.ncbi.nlm.nih.gov/genome)
* ヒトゲノムの主要なデータベース [UCSC Genome Browser](https://genome.ucsc.edu/)
* 主要な生物種を広くカバーするゲノムデータベース [Ensembl](https://www.ensembl.org/), [TogoGenome](http://togogenome.org/)
* ゲノム解読の進行している生物種を網羅したデータベース [JGI GOLD](https://gold.jgi.doe.gov/)

ここでは NCBI Genome を使ってみます。

![NCBI Genome - Takifugu rubripes](images/Genome-Torafugu.png)

学名を入力すると、ゲノムが登録されている場合、DNA 配列やアノテーションおよび関連論文などを参照することができ、さらにゲノムブラウザを利用して遺伝子の分布や構造などの情報を取得することができます。

### 論文データベースPubMedを検索

PubMedは生命医科学論文のほとんどを集積したデータベースです。
これを使って、トラフグのゲノム解読の論文が出ているかどうか調べてみましょう。

![Takifugu rubripes whole genome - PubMed - NCBI](images/PubMed-Torafugu.png)

ゲノムが決まっている場合、学名と genome、whole genome、draft genome、complete genome、genome assembly などのキーワードで論文が見つかってくることが多いです。

PubMed では新しい論文が先頭に表示されますので、ゲノム論文はページをたどって古い論文を見ていかないと見つからないことが多いです。運が良ければ、Best matches の欄にゲノム論文が選ばれて出てくることもあるでしょう。

もし色々なキーワードで探しても見つからない場合は、Google でも検索してみます。
たとえば、スパイシーなゲノム弁当に欲しくなるニンニク Allium sativum L の場合、こんな感じでした。

![Allium sativum genome - Google](images/Google-Garlic.png)

残念ながら、現時点ではまだニンニク本体のゲノム論文は出ていないように見えます。

### Step 5: ゲノム論文とタクソノミー ID を表に追加

もし、新しい食材のゲノムが出ていることが確認できたら、ぜひゲノム弁当の表に追記してください。論文を読んで、ゲノムサイズや遺伝子数、どのチームが解読したかなどが分かればあわせて記載頂けると大変助かります。自信がない場合は DBCLS の片山 <ktym@dbcls.jp> または YCAM の津田さんまでご相談ください。

![GenomeBento 1](images/GenomeBento-1.png)

![GenomeBento 2](images/GenomeBento-2.png)

これで ToyToy さんも安心して新しい食材を加えた、より美味しいお弁当を作ることができるようになります！！

