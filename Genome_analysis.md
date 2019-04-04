# ゲノムのマッピング

ゲノム DNA をシーケンスした、`1_hakusai` と `5_hiyokomame` のリードをそれぞれ白菜とヒヨコ豆のリファレンスゲノム配列にマッピングしてみます。

TODO: 作成時は `1_hakusai_genome.sam` だが、続きの説明が `hakusai.sam` なのを調整

## 白菜の場合

リファレンス配列は `GenomeBento/databases/genomes` フォルダに用意した `3711-hakusai.fasta` を、マッピングするファイルは `MinION/1_hakusai` のフォルダの中の `fastq_pass` にある全ての FASTQ ファイル `*.fastq` をつなげた `samples/1_hakusai.fastq` ファイルを指定しています。

```sh
cd ~/Desktop/GenomeBento
minimap2 -a databases/genomes/3711-hakusai.fasta samples/1_hakusai.fastq > mappings/1_hakusai_genome.sam
```

なお、`minimap2` コマンドが見つからない場合、パスが通ったところにインストールされていない可能性があります。今回は `GenomeBento` フォルダの中にも `minimap2` コマンドを入れてありますので、`./minimap2` のように書き換えると実行できるはずです。

```sh
cd ~/Desktop/GenomeBento
./minimap2 databases/genomes/hakusai samples/1_hakusai.fastq > mappings/1_hakusai_genome.sam
```

## ヒヨコ豆の場合

ヒヨコ豆のデータを使う場合は、白菜の例からデータベース名を `chickpea` に変更し、FASTQ ファイルは `samples/5_hiyokomame.fastq` を使います。

```sh
cd ~/Desktop/GenomeBento
minimap2 -a databases/genomes/chickpea samples/5_hiyokomame.fastq > mapping/5_chickpea_genome.sam
```

ヒヨコ豆の場合、以下の説明では `hakusai` を `chickpea` に読み替えてください。

## IGVを使ってマッピング結果を見てみる

[IGV](http://software.broadinstitute.org/software/igv/) は、米国 Broad Institute で開発されているゲノムブラウザです。リファレンスゲノムおよびマッピングされたリードを表示することができ、ヒトゲノムのように巨大なゲノムにも対応しています。

### IGVを起動 (Dockのアイコンをクリック)

![IGV](images/IGV-icon.png)

はじめて起動した場合、ヒトゲノムのバージョン hg19 が選ばれた状態になっています。

![IGV](images/IGV-initial.png)

### リファレンスゲノム配列を読み込む

他の生物種のゲノムに切り替えるためには、IGV の Genomes メニューから「Load Genome from File...」を選択します。

![IGV](images/IGV-load-genome.png)

白菜だったら `3711-hakusai.genome` など、各グループで読んだ対象の食材の `.genome` ファイルを選択して Open をクリックします。

![Genome file](images/IGV-load-genome-file.png)

読み込めたらこのようにメニューから白菜の各染色体に相当する DNA 配列がメニューに出てきます。

![Loaded Genome](images/IGV-genome-hakusai.png)

ここで読み込んだ `3711-hakusai.genome` のような `.genome` ファイルは[リファレンス配列の準備](Genome_preparation.md)で作成したものです。

### SAM ファイルのソートとインデックス作成

`SAM` ファイルは、できたままだとマッピングされたリードがゲノム座標順に並んでいないので、ゲノムブラウザで見る際に効率が悪くなります。このため、通常はソートしてインデックス作成したものを利用します。後述の `samtools` でもできますが、ここではIGVに内蔵されている `igvtools` を使ってみます。

#### igvtools の起動

Tools メニューから「Run igvtools...」を選んで `igvtools` を起動します。

![igvtools](images/IGV-igvtools.png)

#### igvtools でのソート

Command を `Sort` に変更して、Input File の `Browse` ボタンをクリックして、作成した `hakusai.sam` を選択します。

![igvtools](images/IGV-igvtools-load.png)

`Run` ボタンでソートを実行すると `hakusai.sorted.sam` のように名前に `.sorted` がついた SAM ファイルができます。Done と表示されれば完了です（わりと一瞬で終わるかと思います）。

![igvtools](images/IGV-igvtools-sort.png)

#### igvtools でのインデックス作成

Command を `Index` に変更して、この `hakusai.sorted.sam` ファイルに高速化のためのインデックスを作成します（`Sort` から `Index` に変更しただけだと `hakusai.sorted.sam` ではなく元の `hakusai.sam` ファイルが選ばれている可能性があるので注意）。

![igvtools](images/IGV-igvtools-index.png)

結果、元の `hakusai.sam` ファイルに対して、ソートされた `hakusai.sorted.sam` ファイルと、それに対するインデックスファイル `hakusai.sorted.sam.sai` ができます。

### マッピングの結果を読み込む

「Load from File...」から、作成した `hakusai.sorted.sam` ファイルを選択します。

![IGV](images/IGV-load-sam.png)

あとはズームしていくと、マッピングされたリードが出てきます。

![IGV](images/IGV-result.png)

リードの数が少ない場合、マップされた領域を探すのがちょっと大変ですが、SAM ファイルの中を見ると配列ごとにマッピングされた座標が分かるのでヒントになります。このことから、ゲノムを完全に読むには大量にシーケンスする必要があることが分かると思います。

![SAM](images/IGV-sorted-sam.png)

## BAM ファイルの作成とマップ率の確認

BAM ファイルはSAMファイルをバイナリにして、省スペース化および高速化をはかるためのフォーマットです。SAMファイルをBAMファイルに変換するには [samtools](http://samtools.sourceforge.net/) を使います。`samtools` を使うと、リードのマッピング率なども調べることができます。

### BAM形式に変換

```sh
samtools view -bS hakusai.sam > hakusai.bam
```

### 染色体上の座標でソートする

```sh
samtools sort hakusai.bam hakusai.sorted
```
### インデックスを作成

```sh
samtools index hakusai.sorted.bam
```

### マッピングの状況を見てみる

```sh
samtools flagstat hakusai.sorted.bam
```

今回の白菜の例では下記のような結果になりました。

```
% samtools flagstat hakusai.sorted.bam
9051 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
6067 + 0 mapped (67.03%:nan%)
6067 + 0 paired in sequencing
2999 + 0 read1
3068 + 0 read2
0 + 0 properly paired (0.00%:nan%)
3844 + 0 with itself and mate mapped
2223 + 0 singletons (36.64%:nan%)
3422 + 0 with mate mapped to a different chr
2349 + 0 with mate mapped to a different chr (mapQ>=5)
```

67%のリードがマッピングされたということで、なかなか良かったのではないでしょうか。

TODO: minimap2 でのマップ率に改訂

### bwa コマンドの使い方

TODO: `bwa` から `minimap2` に変更した部分を調整

`bwa` で MinION のロングリードをマッピングするには、サブコマンド `mem` を利用します。コマンドラインオプションは下記のようになります。

```sh
bwa mem リファレンス配列名 入力FASTQファイル > 出力SAMファイル
```

マッピング対象とするリファレンス配列は、あらかじめ準備しておく必要があります。今回は、デスクトップの、`GenomeBento`の中の`databases`フォルダ以下に置いてあります。

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

これらのファイルの作り方に興味のある方は、[リファレンス配列の準備](Genome_preparation.md)を参照してください。
