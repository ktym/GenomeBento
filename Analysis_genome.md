# ゲノムのマッピング

## BWAを用いてリードを参照配列にマッピング

[BWA](http://bio-bwa.sourceforge.net/)は、Burrows-Wheeler Transformation アルゴリズムを用いた高速な配列マッピングを実現しています。Illuminaシーケンサーのショートリードについては様々なマッピングツールがありますが、比較するとMinIONなどのロングリードに対応したツールはまだ少なめです。

`bwa`コマンドは広く使われており十分に高速なのですが、他にも産総研の[LAST](http://last.cbrc.jp/)もロングリードに対応しているほか、東大で開発されている[minialign](https://github.com/ocxtal/minialign)はより高速で正確と謳われており、国内のバイオインフォマティクス研究者の活躍にも期待したいところです。

### BWAコマンドの使い方

`bwa`でMinIONのロングリードをマッピングするには、サブコマンド`bwasw`を利用します。コマンドラインオプションは下記のようになります。

```sh
bwa bwasw リファレンス配列名 入力FASTQファイル > 出力SAMファイル
```

マッピング対象とするリファレンス配列は、あらかじめ準備しておく必要があります。今回は、デスクトップの、GenomeBentoの中のdatabasesフォルダ以下に置いてあります。

```
GenomeBento/databases
|-- genomes
|   |-- banana
|   |-- cabbage
|   |-- carrot
|   |-- chickpea
|   |-- hakusai
|   |-- human
|   |-- rice
|   `-- tomato
`-- rbcL
    |-- rbcL_all
    `-- rbcL_bento
```

実際はそれぞれ、`.amb`, `.ann`, `.bwt`, `.pac`, `.sa` の拡張子がついたファイルの組になります。

#### リファレンス配列の準備

（ここは読み飛ばして構いません）

自分でリファレンス配列のファイルを作る場合は、ゲノム弁当の[食材の追加方法](AdditionalGenomeBento.md)で紹介したように、たとえばNCBI Genomeなどでゲノム配列データのIDを取得します。

![NCBI Genome - Torafugu](images/Genome-Torafugu.png)

ただし、[ココ](https://www.ncbi.nlm.nih.gov/genome/63)↑にリンクされているトラフグの[染色体１番](https://www.ncbi.nlm.nih.gov/nuccore/NC_018890.1)のRefSeqデータベースエントリ`NC_018890.1`には、ゲノム配列自体は含まれていません。

![NCBI RefSeq - Torafugu](images/RefSeq-Torafugu.png)

左上の「FASTA」と書かれたリンクをクリックすると配列データをダウンロードできます。トラフグの全ゲノムを取得するには、この操作を染色体の数繰り返して一つのファイルにまとめます。この作業はちょっと面倒なので、たとえば拙作の[TogoWS](http://togows.org/)を使って工夫してみます。

```sh
for id in NC_018890.1 NC_018891.1 NC_018892.1 NC_018893.1 NC_018894.1 NC_018895.1 NC_018896.1 NC_018897.1 NC_018898.1 NC_018899.1 NC_018900.1 NC_018901.1 NC_018902.1 NC_018903.1 NC_018904.1 NC_018905.1 NC_018906.1 NC_018907.1 NC_018908.1 NC_018909.1 NC_018910.1 NC_018911.1 NC_004299.1
do
  curl http://togows.org/entry/nucleotide/${id}.fasta >> 31033-torafugu.fasta
done
```

少し時間がかかりますが、これで指定したトラフグの各染色体のゲノムDNA配列がすべて`31033-torafugu.fasta`ファイルにFASTA形式でダウンロードできます。

ゲノム配列が入ったFASTAファイルに対し、`bwa index`コマンドを用いてインデックスを作ります。このとき`-p`オプションで`bwa bwasw`にオプションで渡すリファレンス配列名を指定します（FASTAファイル名、リファレンス配列名は自由に決めて構いません）。

```sh
bwa index -p torafugu 31033-torafugu.fasta
```

このようにして作ったのが、今回ゲノム弁当の解読用に用意した上記の`banana`, `cabbage`, `carrot`, `chickpea`, `hakusai`, `rice`, `tomato`などになります。なお、`human`に関しては2019年2月に東北メディカル・メガバンク機構から公開されたばかりの日本人のリファレンスゲノム配列JG1です。実験の途中に間違って自分のDNAが入ってしまって（コンタミして）いないか、確認のために使うことができると思います。

### リファレンス配列へのマッピング

```sh
bwa bwasw databases/genomes/human samples/Shoga/gDNA/FAH64410_caf194bd39f018b128560eba1ade5b8c7f6afc21_* > mappings/human+shoga/human+shoga_gDNA.sam
```

### ここでディレクトリを移動

```sh
cd appings/human+shoga/
```

## BAM形式に変換

```sh
samtools view -bS human+shoga_gDNA.sam > human+shoga_gDNA.bam
```

## 染色体位置でソートする

```sh
samtools sort human+shoga_gDNA.bam human+shoga_gDNA.sorted
```
## Indexを作成

```sh
samtools index human+shoga_gDNA.sorted.bam
```

## マッピングの状況を見てみる

```sh
samtools flagstat human+shoga_gDNA.sorted.bam
```
