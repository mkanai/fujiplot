FROM alpine

ENV version 0.69-6-kanai

# ADD http://circos.ca/distribution/circos-${version}.tgz /tmp/
COPY circos-0.69-6-kanai.tgz /tmp/

RUN     cd /opt/ \
        && tar xzvf /tmp/circos-0.69-6-kanai.tgz \
        && mv circos-${version} circos \
        && sed -i 's/max_points_per_track = 25000/max_points_per_track = 100000/' /opt/circos/etc/housekeeping.conf \
        && apk add --update --no-cache perl gd jpeg freetype \
        && apk add --update --no-cache --virtual=deps make gd-dev jpeg-dev freetype-dev apkbuild-cpan gcc  musl-dev perl-dev \
        && wget -O - http://cpanmin.us | perl - --self-upgrade  \
        && cpanm Math::Bezier Math::Round Readonly::Tiny Readonly Config::General Params::Validate Font::TTF::Font Regexp::Common Math::VecStat Text::Format SVG Clone List::MoreUtils  \
        && cpanm -force GD Number::Format \
        && cpanm Statistics::Basic Set::IntSpan \
        && cpanm -force Try::Tiny \
        && rm -rf /var/cache/apk/* /tmp/circos-${version}.tgz \
        && apk del deps


CMD ["-help"]
ENTRYPOINT ["/opt/circos/bin/circos"]
