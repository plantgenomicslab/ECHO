
#!/usr/bin/env python
import sys
import argparse
import re

parser = argparse.ArgumentParser(description='Usage:')

parser.add_argument('augustus_gff3', help='augustus.gff3')
parser.add_argument('transfrag_gff3', help='transfrag.gff3')
parser.add_argument('genewise_gff3', help='genewise.gff3')
parser.add_argument('intron_gff', help='intron.gff')
parser.add_argument('-o', '--overlap', type=int, default=30, help='overlap')
parser.add_argument('-m', '--min_augustus_transcriptSupport_percentage', type=float, default=10.0, help='min_augustus_transcriptSupport_percentage')
parser.add_argument('-n', '--min_augustus_intronSupport_number', type=int, default=1, help='min_augustus_intronSupport_number')
parser.add_argument('-r', '--min_augustus_intronSupport_ratio', type=float, default=0.01, help='min_augustus_intronSupport_ratio')
parser.add_argument('--more_strict', action='store_true', help='more_strict')

args = parser.parse_args()

augustus = {}
transfrag = {}
genewise = {}
intron = {}
geneRegion = {}
aug_gene = {}
aug_gene_list = []

file_list = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]]

for i,file in enumerate(file_list):
    with open(file) as f:
        for line in f:
            if line.startswith("#") or line.isspace():
                continue
            if "gene" in line and "ID=" in line:
                id = re.search("ID=([^;\s]+)", line).group(1)
                _ = line.strip().split("\t")
                geneRegion[f"{_[3]}\t{_[4]}"][id] = 1
            if i==0:
                if "mRNA" in line and "ID=" in line:
                    id = re.search("ID=([^;\s]+)", line).group(1)
                    aug_gene[aug_gene][id] = 1
                    _ = line.strip().split("\t")
                    geneRegion[f"{_[3]}\t{_[4]}"][id] = 1
                augustus[id] += line
            elif i==1:
                transfrag[id] += line
            elif i==2:
                genewise[id] += line
            elif i==3:
                if "intron" in line:
                    _ = line.strip().split("\t")
                    intron[f"{_[0]}\t{_[3]}\t{_[4]}"] = 1

gene_out = {}
for region in sorted(region.keys()):
    ID = list(region[region].keys())
    ID_str = ','.join(ID)

    status_augustus, status_transfrag, status_genewise = (0, 0, 0)
    for id in ID:
        if id in augustus:
            status_augustus = 1
        if id in transfrag:
            status_transfrag = 1
        if id in genewise:
            status_genewise = 1

    if status_augustus == 1:
        for gene_id in ID:
            if gene_id in augustus:
                intron_support = intron_support(gene_id)
                gff3_out = add_UTR(gene_id, ID)

                for cds in gff3_out[gene_id]["CDS"]:
                    gene_out[gene_id]["CDS"][cds] = 1
                for exon in gff3_out[gene_id]["exon"]:
                    gene_out[gene_id]["exon"][exon] = 1

                gene_out[gene_id]["strand"] = gff3_out[gene_id]["strand"]
                gene_out[gene_id]["chr"] = gff3_out[gene_id]["chr"]
                gene_out[gene_id]["intron_support"] = intron_support

                if intron_support.startswith("0"):
                    gene_out[gene_id]["pfam"] = 1
    elif status_transfrag == 1:
        for gene_id in ID:
            if gene_id in transfrag:
                gff3_out = combine_transfrag_and_genewise(gene_id, ID)
                for cds in gff3_out[gene_id]["CDS"]:
                    gene_out[gene_id]["CDS"][cds] = 1
                for exon in gff3_out[gene_id]["exon"]:
                    gene_out[gene_id]["exon"][exon] = 1

                gene_out[gene_id]["strand"] = gff3_out[gene_id]["strand"]
                gene_out[gene_id]["chr"] = gff3_out[gene_id]["chr"]
    elif status_genewise == 1:
        for gene_id in ID:
            if gene_id in genewise:
                genewise_info = genewise[gene_id]
                for line in genewise_info.split("\n"):
                    fields = line.split("\t")
                    if fields[2] == "CDS":
                        gene_out[gene_id]["CDS"][f"{fields[3]}\t{fields[4]}\t{fields[5]}\t{fields[6]}\t{fields[7]}"] = 1
                        gene_out[gene_id]["exon"][f"{fields[3]}\t{fields[4]}"] = 1
                    elif fields[2] == "gene":
                        gene_out[gene_id]["chr"] = fields[0]
                        gene_out[gene_id]["strand"] = fields[6]


for gene_id in aug_gene:
    out = {}
    gene_pos = {}
    if_intron_support = None

    fields = aug_gene_info[gene_id].split("\t")
    chr, strand, gene_score = fields[0], fields[6], fields[5]
    augustus_transcriptSupport_percentage, augustus_intronSupport_number, augustus_intronSupport_ratio, augustus_intronSupport = None, None, None, None
    if re.search(r"exonHintRatio=([^;\s]+)", fields[8]):
        augustus_transcriptSupport_percentage = re.search(r"exonHintRatio=([^;\s]+)", fields[8]).group(1)
    if re.search(r"intronSupport=(\d+)\/(\d+)", fields[8]):
        match = re.search(r"intronSupport=(\d+)\/(\d+)", fields[8])
        augustus_intronSupport_number = match.group(1)
        augustus_intronSupport = match.group(0)
        if match.group(2) == 0:
            augustus_intronSupport_ratio = 0
        else:
            augustus_intronSupport_ratio = match.group(1) / match.group(2)

    mRNA_number = 0
    for mRNA_id in sorted(aug_gene[gene_id].keys()):
        if_intron_support = 1 if mRNA_id not in gene_out or "pfam" not in gene_out[mRNA_id]
        out = {}
        position = {}
        mRNA_number += 1

        if strand == "+":
            number = 0
            for key in sorted(gene_out[mRNA_id]["CDS"].keys()):
                number += 1
                fields = key.split("\t")
                position[fields[0]] = 1
                gene_pos[fields[0]] = 1
                position[fields[1]] = 1
                gene_pos[fields[1]] = 1
                out[f"{chr}\t.\tCDS\t{key}\tID={mRNA_id}.CDS{number};Parent={mRNA_id};"] = fields[0]

            number = 0
            for key in sorted(gene_out[mRNA_id]["exon"].keys(), reverse=True):
                number += 1
                fields = key.split("\t")
                position[fields[0]] = 1
                gene_pos[fields[0]] = 1
                position[fields[1]] = 1
                gene_pos[fields[1]] = 1
                out[f"{chr}\t.\texon\t{key}\t.\t{strand}\t.\tID={mRNA_id}.exon{number};Parent={mRNA_id};"] = fields[0]


def intron_support(gene_id):
    augustus_info = augustus[gene_id]
    augustus_info = augustus_info.split("\n")
    total_num, support_num = 0, 0
    for line in augustus_info:
        if "\tintron\t" in line:
            total_num += 1
            data = line.split("\t")
            if f"{data[0]}\t{data[3]}\t{data[4]}" in intron:
                support_num += 1
    return f"{support_num}/{total_num}"



def add_UTR(augustus_id, *id):
    augustus_info = augustus[augustus_id]
    gff3_out, gene_id, strand, augustus_CDS = {}, None, None, []

    match = re.search(r"(\S+?)\t\S+?\tmRNA\t\d+?\t\d+?\t\S+?\t(\S+?)\t\S+?\tID=([^;\s]+)", augustus_info)
    if match:
        gene_id = match.group(3)
        strand = match.group(2)
        gff3_out[gene_id] = {"chr": match.group(1), "strand": strand}

    for line in augustus_info.split("\n"):
        data = line.split("\t")
        if data[2] == "start_codon":
            gff3_out[gene_id]["start_codon"] = {f"{data[3]}\t{data[4]}\t{data[5]}\t{data[6]}\t{data[7]}": 1}
        elif data[2] == "stop_codon":
            gff3_out[gene_id]["stop_codon"] = {f"{data[3]}\t{data[4]}\t{data[5]}\t{data[6]}\t{data[7]}": 1}
        elif data[2] == "CDS":
            gff3_out[gene_id]["CDS"] = {f"{data[3]}\t{data[4]}\t{data[5]}\t{data[6]}\t{data[7]}": 1}
            augustus_CDS.append(f"{data[3]}\t{data[4]}")
    augustus_CDS.sort()

    utr5_cds = augustus_CDS[0].split("\t")
    utr3_cds = augustus_CDS[-1].split("\t")

    augustus_CDS.pop(0)
    augustus_CDS.pop(-1)
    for exon in augustus_CDS:
        data = exon.split("\t")
        gff3_out[gene_id]["exon"] = {f"{data[0]}\t{data[1]}": 1}

    utr5_status, utr3_status = 0, 0
    for transfrag_id in id:
        if transfrag_id in transfrag:
            transfrag_info = transfrag[transfrag_id]
            utr5_ok, utr3_ok = 0, 0
            transfrag_exon = []
            for line in transfrag_info.split("\n"):
                data = line.split("\t")
                if data[2] == "CDS":
                    if int(data[3]) <= int(utr5_cds[0]) and int(data[4]) >= int(utr5_cds[1]):
                        utr5_ok = 1
                    if int(data[3]) <= int(utr3_cds[0]) and int(data[4]) >= int(utr3_cds[1]):
                        utr3_ok = 1
                elif data[2] == "exon":
                    transfrag_exon.append(f"{data[3]}\t{data[4]}")
            
            transfrag_exon.sort()
            
            if utr5_ok == 1:
                for exon in transfrag_exon:
                    data = exon.split("\t")
                    if int(data[0]) <= int(utr5_cds[1]):
                        if int(data[1]) > int(utr5_cds[1]):
                            gff3_out[gene_id]["exon"] = {f"{data[0]}\t{utr5_cds[1]}": 1}
                        else:
                            gff3_out[gene_id]["exon"] = {f"{data[0]}\t{data[1]}": 1}
                utr5_status = 1
            
            if utr3_ok == 1:
                for exon in transfrag_exon:
                    data = exon.split("\t")
                    if int(data[1]) >= int(utr3_cds[0]):
                        if int(data[0]) < int(utr5_cds[0]):
                            gff3_out[gene_id]["exon"] = {f"{utr5_cds[0]}\t{data[1]}": 1}
                        else:
                            gff3_out[gene_id]["exon"] = {f"{data[0]}\t{data[1]}": 1}
                utr3_status = 1
    
    if utr5_status == 0:
        gff3_out[gene_id]["exon"] = {f"{utr5_cds[0]}\t{utr5_cds[1]}": 1}
    if utr3_status== 0:
        gff3_out[gene_id]["exon"] = {f"{utr3_cds[0]}\t{utr3_cds[1]}": 1}

    print(f"{gff3_out[gene_id]['chr']}\t.\tstart_codon\t{gff3_out[gene_id]['start_codon']}\tID={gene_id}.t1.start_codon;Parent={gene_id}.t1")
    print(f"{gff3_out[gene_id]['chr']}\t.\tstop_codon\t{gff3_out[gene_id]['stop_codon']}\tID={gene_id}.t1.stop_codon;Parent={gene_id}.t1")
    for cds in sorted(gff3_out[gene_id]["CDS"]):
        print(f"{gff3_out[gene_id]['chr']}\t.\tCDS\t{cds}\tID={gene_id}.t1.cds;Parent={gene_id}.t1")
    for exon in sorted(gff3_out[gene_id]["exon"]):
        print(f"{gff3_out[gene_id]['chr']}\t.\texon\t{exon}\t.\t{gff3_out[gene_id]['strand']}\t.\tID={gene_id}.t1.exon;Parent={gene_id}.t1")

    return gff3_out

#############
foreach my $gene_id (sort keys %gene_out) {
    my %out;
    my $chr = $gene_out{$gene_id}{'chr'};
    my $strand = $gene_out{$gene_id}{'strand'};
    my %position;

    if ($strand eq "+") {
        my $number = 0;
        foreach (sort {$a <=> $b} keys %{$gene_out{$gene_id}{"CDS"}}) {
            $number ++;
            @_ = split /\t/; 
            $position{$_[0]} = 1;
            $position{$_[1]} = 1;
            $out{"$chr\t\.\tCDS\t$_\tID=$gene_id.t1.CDS$number;Parent=$gene_id.t1;"} = $_[0];
        }
        $number = 0;
        foreach (sort {$a <=> $b} keys %{$gene_out{$gene_id}{"exon"}}) {
            $number ++;
            @_ = split /\t/;
            $position{$_[0]} = 1;
            $position{$_[1]} = 1;
            $out{"$chr\t\.\texon\t$_\t\.\t$strand\t\.\tID=$gene_id.t1.exon$number;Parent=$gene_id.t1;"} = $_[0];
        }
    }
    elsif ($strand eq "-") {
        my $number = 0;
        foreach (sort {$b <=> $a} keys %{$gene_out{$gene_id}{"CDS"}}) {
            $number ++;
            @_ = split /\t/; 
            $position{$_[0]} = 1;
            $position{$_[1]} = 1;
            $out{"$chr\t\.\tCDS\t$_\tID=$gene_id.t1.CDS$number;Parent=$gene_id.t1;"} = $_[0];
        }
        $number = 0;
        foreach (sort {$b <=> $a} keys %{$gene_out{$gene_id}{"exon"}}) {
            $number ++;
            @_ = split /\t/;
            $position{$_[0]} = 1;
            $position{$_[1]} = 1;
            $out{"$chr\t\.\texon\t$_\t\.\t$strand\t\.\tID=$gene_id.t1.exon$number;Parent=$gene_id.t1;"} = $_[0];
        }
    }

    foreach (keys %position) { delete $position{$_} unless $_; }
    my @position = sort {$a <=> $b} keys %position;
    print "$chr\t\.\tgene\t$position[0]\t$position[-1]\t\.\t$strand\t\.\tID=$gene_id;\n";
    print "$chr\t\.\tmRNA\t$position[0]\t$position[-1]\t\.\t$strand\t\.\tID=$gene_id.t1;Parent=$gene_id;\n";

    my @out;
    if ($strand eq "+") {
        @out = sort {$out{$a} <=> $out{$b} or $b cmp $a} keys %out;
    }
    elsif ($strand eq "-") {
        @out = sort {$out{$b} <=> $out{$a} or $b cmp $a} keys %out;
    }
    foreach (@out) {
        print "$_\n";
    }
}

sub combine_transfrag_and_genewise {
    # 输入的第一个参数是transfrag的基因ID
    # 输入的第二个参数是数组：该基因区域的所有基因ID
    my $gene_id = shift @_;
    my $transfrag_info = $transfrag{$gene_id};
    my @id = @_;

    # 得到chr，strand和transfrag基因的完整性信息
    my (%gff3_out, $strand, $intergrity);
    if ($transfrag_info =~ m/(\S+?)\t\S+?\tgene\t\d+?\t\d+?\t\S+?\t(\S+?).*?Integrity=([^;\s]+)/) {
        $strand = $2;
        $gff3_out{$gene_id}{"chr"} = $1;
        $gff3_out{$gene_id}{"strand"} = $strand;
        $intergrity=$3;
    }

    # 得到transfrag基因的CDS和exon信息
    my (@CDS, %CDS, %exon);
    foreach my $line (split /\n/, $transfrag_info) {
        @_ = split /\t/, $line;
        if ($_[2] eq "CDS") {
            push @CDS, "$_[3]\t$_[4]\t$_[7]";
            $CDS{"$_[3]\t$_[4]\t$_[7]"} = 1;
        }
        elsif ($_[2] eq "exon") {
            $exon{"$_[3]\t$_[4]"} = 1;
        }
    }
    @CDS = sort {$a <=> $b} @CDS;

    # transfrag基因是5prime_partial,3prime_partial,internal,strand等情况下，分别使用genewise结果进行整合
    # transfrag基因首尾CDS和genewise的CDS有重叠，且具有相同的frame，则进行整合。
    # transfrag基因是complete则不整合
    if ($intergrity == '5prime_partial' or $intergrity == 'internal') {
        if ($strand eq "+") {
            my @cds_first = split /\t/, $CDS[0];
            foreach (@id) {
                if (exists $genewise{$_}) {
                    my $genewise_info = $genewise{$_};
                    my @genewise_cds;
                    foreach (split /\n/, $genewise_info) {
                        @_ = split /\t/,
                        push @genewise_cds, "$_[3]\t$_[4]\t$_[7]" if $_[2] eq 'CDS';
                    }
                    @genewise_cds = sort {$a <=> $b} @genewise_cds;

                    my (@cds_ok, $merge);
                    foreach (@genewise_cds) {
                        @_ = split /\t/;
                        if ($_[0] < $cds_first[0] && $_[1] >= $cds_first[0]) {
                            my $frame = ($cds_first[0] - $_[0] - $_[2]) % 3;
                            $frame = 1 if $frame == 2;
                            $frame = 2 if $frame == 1;
=cut
       =================        cds_first
   21021021
   12345678
     ^^^###
   =================           genewise
   2

0  =>  0
1  =>  2
2  =>  1

       ==================   cds_first
   1021021
   12345678
    ^^^###
   =================           genewise
=cut

                            if ($frame == $cds_first[2]) {
                                push @cds_ok, "$_[0]\t$cds_first[1]\t$_[2]";
                                $merge = 1;
                                last;
                            }
                        }
                        else {
                            push @cds_ok, $_;
                        }
                    }

                    if ($merge == 1) {
                        delete $CDS{$CDS[1]};
                        foreach (@cds_ok) {
                            $CDS{$_} = 1;
                        }
                    }
                }
            }
        }
        elsif ($strand eq "-") {
            my @cds_first = split /\t/, $CDS[-1];
            foreach (@id) {
                if (exists $genewise{$_}) {
                    my $genewise_info = $genewise{$_};
                    my @genewise_cds;
                    foreach (split /\n/, $genewise_info) {
                        @_ = split /\t/,
                        push @genewise_cds, "$_[3]\t$_[4]\t$_[7]" if $_[2] eq 'CDS';
                    }
                    @genewise_cds = sort {$b <=> $a} @genewise_cds;
                    my (@cds_ok, $merge);
                    foreach (@genewise_cds) {
                        @_ = split /\t/;
                        if ($_[1] > $cds_first[1] && $_[0] <= $cds_first[1]) {
                            my $frame = ($_[1] - $cds_first[1] - $_[2]) % 3;
                            $frame = 2 if $frame == 1;
                            $frame = 1 if $frame == 2;
=cut
          =================                  cds_fisrt
                       12012012
                       12345678
                       ###^^^
              =================              genewise
                              2
=cut

                            if ($frame == $cds_first[2]) {
                                push @cds_ok, "$cds_first[0]\t$_[1]\t$_[2]";
                                $merge = 1;
                                last;
                            }
                        }
                        else {
                            push @cds_ok, $_;
                        }
                    }

                    if ($merge == 1) {
                        delete $CDS{$CDS[-1]};
                        foreach (@cds_ok) {
                            $CDS{$_} = 1;
                        }
                    }
                }
            }
        }
    }
    
    if ($intergrity == '3prime_partial' or $intergrity == 'internal') {
        if ($strand eq "+") {
            my @cds_last = split /\t/, $CDS[-1];
            foreach (@id) {
                if (exists $genewise{$_}) {
                    my $genewise_info = $genewise{$_};
                    my @genewise_cds;
                    foreach (split /\n/, $genewise_info) {
                        @_ = split /\t/,
                        push @genewise_cds, "$_[3]\t$_[4]\t$_[7]" if $_[2] eq 'CDS';
                    }
                    @genewise_cds = sort {$b <=> $a} @genewise_cds;
                    my (@cds_ok, $merge);
                    foreach (@genewise_cds) {
                        @_ = split /\t/;
                        if ($_[1] > $cds_last[1] && $_[0] <= $cds_last[1]) {
                            if ($_[0] > $cds_last[0]) {
                                my $frame = ($_[0] - $cds_last[0] - $cds_last[2]) % 3;
                                $frame = 1 if $frame == 2;
                                $frame = 2 if $frame == 1;
=cut
       2
       ================>               cds_last
       12345678
       21021021
         ^^^###
           ==================>         genewise
=cut
                                if ($frame == $_[2]) {
                                    push @cds_ok, "$cds_last[0]\t$_[1]\t$cds_last[2]";
                                    $merge = 1;
                                    last;
                                }
                            }
                            else {
                                my $frame = ($cds_last[0] - $_[0] - $_[2]) % 3;
                                $frame = 1 if $frame == 2;
                                $frame = 2 if $frame == 1;
                                if ($frame == $cds_last[2]) {
                                    push @cds_ok, "$cds_last[0]\t$_[1]\t$cds_last[2]";
                                    $merge = 1;
                                    last;
                                }
                            }
                        }
                        else {
                            push @cds_ok, $_;
                        }
                    }

                    if ($merge == 1) {
                        delete $CDS{$CDS[-1]};
                        foreach (@cds_ok) {
                            $CDS{$_} = 1;
                        }
                    }
                }
            }
        }
        elsif ($strand eq "-") {
            my @cds_last = split /\t/, $CDS[0];
            foreach (@id) {
                if (exists $genewise{$_}) {
                    my $genewise_info = $genewise{$_};
                    my @genewise_cds;
                    foreach (split /\n/, $genewise_info) {
                        @_ = split /\t/,
                        push @genewise_cds, "$_[3]\t$_[4]\t$_[7]" if $_[2] eq 'CDS';
                    }
                    @genewise_cds = sort {$a <=> $b} @genewise_cds;
                    my (@cds_ok, $merge);
                    foreach (@genewise_cds) {
                        @_ = split /\t/;
                        if ($_[0] < $cds_last[0] && $_[1] >= $cds_last[0]) {
                            if ($_[1] > $cds_last[1]) {
                                my $frame = ($_[1] - $cds_last[1] - $_[2]) % 3;
                                $frame = 1 if $frame == 2;
                                $frame = 2 if $frame == 1;
=cut
            <============             cds_last
                     12012012
                     12345678
       <=====================              genewise
                            2
=cut
                                if ($frame == $cds_last[2]) {
                                    push @cds_ok, "$_[0]\t$cds_last[1]\t$cds_last[2]";
                                    $merge = 1;
                                    last;
                                }
                            }
                            else {
                                my $frame = ($cds_last[1] - $_[1] - $cds_last[2]) % 3;
                                $frame = 1 if $frame == 2;
                                $frame = 2 if $frame == 1;
                                if ($frame == $_[2]) {
                                    push @cds_ok, "$_[0]\t$cds_last[1]\t$cds_last[2]";
                                    $merge = 1;
                                    last;
                                }
                            }
                        }
                        else {
                            push @cds_ok, $_;
                        }
                    }

                    if ($merge == 1) {
                        delete $CDS{$CDS[0]};
                        foreach (@cds_ok) {
                            $CDS{$_} = 1;
                        }
                    }
                }
            }
        }
    }

    if ($intergrity == 'complete') {
        foreach (@CDS) {
            @_ = split /\t/;
            $gff3_out{$gene_id}{"CDS"}{"$_[0]\t$_[1]\t\.\t$strand\t$_[2]"} = 1;
        }
        foreach (keys %exon) {
            $gff3_out{$gene_id}{"exon"}{$_} = 1;
        }
    }
    else {
        my %cds_exon;
        foreach (@CDS) {
            @_ = split /\t/;
            $gff3_out{$gene_id}{"CDS"}{"$_[0]\t$_[1]\t\.\t$strand\t$_[2]"} = 1;
            $cds_exon{"$_[0]\t$_[1]"} = 1;
        }
        foreach (keys %exon) {
            my ($start, $end) = split /\t/;
            my $keep = 1;
            foreach (keys %cds_exon) {
                @_ = split /\t/;
                if ($_[0] <= $end && $_[1] >= $start) {
                    $keep = 0;
                    if ($_[0] >= $start && $_[1] <= $end) {
                        delete $cds_exon{$_};
                        $cds_exon{"$start\t$end"} = 1;
                    }
                }
            } 
            $cds_exon{"$start\t$end"} = 1 if $keep == 1;
        }
        foreach (keys %cds_exon) {
            $gff3_out{$gene_id}{"exon"}{$_} = 1;
        }
    }

    return %gff3_out;
}

sub add_UTR {
    my $augustus_id = shift @_;
    my @id = @_;

    my $augustus_info = $augustus{$augustus_id};
    my (%gff3_out, $gene_id, $strand, @augustus_CDS);

    if ($augustus_info =~ m/(\S+?)\t\S+?\tmRNA\t\d+?\t\d+?\t\S+?\t(\S+?)\t\S+?\tID=([^;\s]+)/) {
        $gene_id = $3;
        $strand = $2;
        $gff3_out{$gene_id}{"chr"} = $1;
        $gff3_out{$gene_id}{"strand"} = $strand;
        #print "$1\t$gene_id\t$strand\n";
    }

    foreach my $line (split /\n/, $augustus_info) {
        @_ = split /\t/, $line;
        if ($_[2] eq "start_codon") {
            $gff3_out{$gene_id}{"start_codon"}{"$_[3]\t$_[4]\t$_[5]\t$_[6]\t$_[7]"} = 1;
        }
        elsif ($_[2] eq "stop_codon") {
            $gff3_out{$gene_id}{"stop_codon"}{"$_[3]\t$_[4]\t$_[5]\t$_[6]\t$_[7]"} = 1;
        }
        elsif ($_[2] eq "CDS") {
            $gff3_out{$gene_id}{"CDS"}{"$_[3]\t$_[4]\t$_[5]\t$_[6]\t$_[7]"} = 1;
            push @augustus_CDS, "$_[3]\t$_[4]";
        }
    }

    @augustus_CDS = sort {$a <=> $b} @augustus_CDS;

    my @utr5_cds = split /\t/, $augustus_CDS[0];
    my @utr3_cds = split /\t/, $augustus_CDS[-1];

    shift @augustus_CDS;
    pop @augustus_CDS;
    foreach (@augustus_CDS) {
        @_ = split /\t/;
        $gff3_out{$gene_id}{"exon"}{"$_[0]\t$_[1]"} = 1;
    }

    my ($utr5_status, $utr3_status) = (0, 0);
    foreach (@id) {
        if (exists $transfrag{$_}) {
            my $transfrag_info = $transfrag{$_};
            my ($utr5_ok, $utr3_ok) = (0, 0);
            my @transfrag_exon;
            foreach my $line (split /\n/, $transfrag_info) {
                @_ = split /\t/, $line;
                if ($_[2] eq "CDS") {
                    if ($_[3] <= $utr5_cds[0] && $_[4] >= $utr5_cds[1]) {
                        $utr5_ok = 1;
                    }
                    if ($_[3] <= $utr3_cds[0] && $_[4] >= $utr3_cds[1]) {
                        $utr3_ok = 1;
                    }
                }
                elsif ($_[2] eq "exon") {
                    push @transfrag_exon, "$_[3]\t$_[4]";
                }
            }

            @transfrag_exon = sort {$a <=> $b} @transfrag_exon;

            if ($utr5_ok == 1) {
                foreach (@transfrag_exon) {
                    @_ = split /\t/;
                    if ($_[0] <= $utr5_cds[1]) {
                        if ($_[1] > $utr5_cds[1]) {
                            $gff3_out{$gene_id}{"exon"}{"$_[0]\t$utr5_cds[1]"} = 1;
                        }
                        else {
                            $gff3_out{$gene_id}{"exon"}{"$_[0]\t$_[1]"} = 1;
                        }
                    }
                }
                $utr5_status = 1;
            }

            if ($utr3_ok == 1) {
                foreach (@transfrag_exon) {
                    @_ = split /\t/;
                    if ($_[1] >= $utr3_cds[0]) {
                        if ($_[0] < $utr5_cds[0]) {
                            $gff3_out{$gene_id}{"exon"}{"$utr5_cds[0]\t$_[1]"} = 1;
                        }
                        else {
                            $gff3_out{$gene_id}{"exon"}{"$_[0]\t$_[1]"} = 1;
                        }
                    }
                }
                $utr3_status = 1;
            }
        }
    }

    if ($utr5_status == 0) {
        $gff3_out{$gene_id}{"exon"}{"$utr5_cds[0]\t$utr5_cds[1]"} = 1;
    }
    if ($utr3_status== 0) {
        $gff3_out{$gene_id}{"exon"}{"$utr3_cds[0]\t$utr3_cds[1]"} = 1;
    }
=cut
    print "$gff3_out{$gene_id}{'chr'}\t\.\tstart_codon\t$gff3_out{$gene_id}{'start_codon'}\tID=$gene_id.t1.start_codon;Parent=$gene_id.t1\n";
    print "$gff3_out{$gene_id}{'chr'}\t\.\tstop_codon\t$gff3_out{$gene_id}{'stop_codon'}\tID=$gene_id.t1.stop_codon;Parent=$gene_id.t1\n";
    foreach (sort {$a <=> $b} keys %{$gff3_out{$gene_id}{"CDS"}}) {
        print "$gff3_out{$gene_id}{'chr'}\t\.\tCDS\t$_\tID=$gene_id.t1.cds;Parent=$gene_id.t1\n";
    }
    foreach (sort {$a <=> $b} keys %{$gff3_out{$gene_id}{"exon"}}) {
        print "$gff3_out{$gene_id}{'chr'}\t\.\texon\t$_\t\.\t$gff3_out{$gene_id}{'strand'}\t\.\tID=$gene_id.t1.exon;Parent=$gene_id.t1\n";
    }
=cut
    return %gff3_out;
}

sub intron_support {
    my $gene_id = $_[0];
    my $augustus_info = $augustus{$gene_id};
    my @augustus_info = split /\n/, $augustus_info;
    my ($total_num, $support_num) = (0, 0);
    foreach (@augustus_info) {
        if (m/\tintron\t/) {
            $total_num ++;
            @_ = split /\t/;
            if (exists $intron{"$_[0]\t$_[3]\t$_[4]"}) {
                $support_num ++;
            }
        }
    }
    #print STDERR "$gene_id\t$support_num\t$total_num\n";
    return "$support_num/$total_num";
}
