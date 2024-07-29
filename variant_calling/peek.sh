for f in $(  find PE/05_filter/ -name "*bcf" );do
    vt peek ${f} -y ${f}.peek.pdf &> ${f}.peek
done
