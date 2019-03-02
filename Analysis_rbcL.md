# rbcL によるメタ植物ゲノム解析

料理については複数の食材が入っているため、野菜のゲノムDNA中から rbcL 遺伝子の領域だけを PCR で増幅して、含まれている食材が検出できるか試してみます。この rbcL は植物が二酸化炭素から炭素を固定するルビスコ (rubisco) という大変重要な酵素（のlarge subunit）を構成する遺伝子で、正確には「リブロース1,5-ビスリン酸カルボキシラーゼ/オキシゲナーゼ」英語では「ribulose 1,5-bisphosphate carboxylase/oxygenase (を略してRuBisCO)」というそうです（リブロースといっても牛肉ではありません…）。なお rbcL は核ゲノムではなく、葉緑体ゲノムに含まれる遺伝子で、全ての植物から見つかるため、DNAによる植物種検出のバーコードとして利用されています。ルビスコの立体構造についての解説や、酵素の働きが鈍いため葉緑体には大量のルビスコが含まれており、その結果世界で一番豊富に存在する酵素となっている、といった興味深い話は [PDBj 今月の分子 011](https://pdbj.org/mom/11) で読むことができます。


- `databases/rbcL` に必要なファイルがあります
![filelist](https://github.com/ktym/GenomeBento/blob/master/images/imag.png "サンプル")

- rbcLライブラリに、リードをマッピング
```
bwa bwasw ~/Desktop/GenomeBento/databases/rbcL_all ../FAH64471_43586857e208cde4192344268331332d0af3dbf0_0.fastq 
```
- マッピング結果を見てみる（各リードがどのrbcLに張り付いたのかを確認）
```
./check.sh mappings/rbcL+shoga_rbcL/rbcL+shoga.sam | less
./check.sh mappings/rbcL+shoga_rbcL/rbcL+shoga.sam | grep YCAMGB_
```

