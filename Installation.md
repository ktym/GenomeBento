# MinION のリードを参照配列にマッピング

Illumina シーケンサーのショートリードについては様々なマッピングツールがありますが、MinION などのロングリードに対応したツールはまだ少なめです。[BWA](http://bio-bwa.sourceforge.net/) は、Burrows-Wheeler Transformation アルゴリズムを用いた高速な配列マッピングを実現しており、ショートリードのアラインメントで長く使われてきました。`bwa` にはロングリードをアラインメントするための `bwa bwasw` や改良版 `bwa mem` などのサブコマンドがあります。日本国内の研究者も活躍しており、産総研の [LAST](http://last.cbrc.jp/) もロングリードに対応しているほか、東大で開発されている [minialign](https://github.com/ocxtal/minialign) はより高速で正確と謳われています。現在のところ `minialign` の作者の鈴木さんによると [minimap2](https://github.com/lh3/minimap2) がベストな選択だそうです。

マッピングしたデータは SAM 形式のファイルになります。これを扱うためのツールには [samtools](http://samtools.sourceforge.net/) があります。また、結果を可視化するためには [IGV](https://software.broadinstitute.org/software/igv/) などのブラウザを利用します。これらのツールをインストールしておきましょう。

## `minimap2` コマンドのインストール

`minimap2` をインストールするには、[リリース](https://github.com/lh3/minimap2/releases)のページから最新版を取得します。バージョン番号は適宜最新のものに読み替えてください。いずれの場合も、動作が確認できたら `minimap2` コマンドをパスの通ったディレクトリ（`/usr/local/bin` など）にコピーします。

### Linux の場合

Linux環境であれば `minimap2-2.16_x64-linux.tar.bz2` をダウンロードして解凍するとすぐに利用できる実行ファイルが入っています。

```sh
curl -O https://github.com/lh3/minimap2/releases/download/v2.16/minimap2-2.16_x64-linux.tar.bz2
tar xvf minimap2-2.16_x64-linux.tar.bz2
./minimap2-2.16_x64-linux/minimap2
Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]
  :
```

たとえば /usr/local/bin/minimap2 としてインストールする場合は下記のようになります。

```sh
sudo cp ./minimap2-2.16_x64-linux/minimap2 /usr/local/bin
```

### Mac の場合

Mac OS X の場合はソースコードを取得してコンパイルする必要があります。

```sh
curl -O https://github.com/lh3/minimap2/releases/download/v2.16/minimap2-2.16.tar.bz2
tar xvf minimap2-2.16.tar.bz2
cd minimap2-2.16
make
./minimap2
Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]
  :
```

たとえば /usr/local/bin/minimap2 としてインストールする場合は下記のようになります。

```sh
sudo cp ./minimap2 /usr/local/bin
```

### Windows の場合

Windows 環境でのインストール方法は環境構築が必要になるためここでは紹介しません（スミマセン）。

## `bwa` コマンドのインストール

YCAM の実習では `minimap2` の代わりに `bwa` を用いました。このため `bwa` コマンドのインストール方法も記載しておきます。

TODO: 追記する

## `samtools` コマンドのインストール

TODO: 追記する

## `IGV` のインストール

TODO: 追記する
