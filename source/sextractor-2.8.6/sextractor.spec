%define name sextractor
%define version 2.8.6
%define release 1
%undefine _missing_build_ids_terminate_build

Summary: Extract catalogs of sources from astronomical images
Name: %{name}
Version: %{version}
Release: %{release}
Source0: ftp://ftp.iap.fr/pub/from_users/bertin/sextractor/%{name}-%{version}.tar.gz
URL: http://astromatic.iap.fr/software/%{name}
License: CeCILL
Group: Sciences/Astronomy
BuildRoot: %{_tmppath}/%{name}-buildroot
BuildRequires: pkgconfig
BuildRequires: fftw-devel >= 3.1
BuildRequires: atlas-devel >= 3.6.0

%description
SExtractor stands for ``Source Extractor'': a software for making catalog of sources from astronomical images.

%prep
%setup -q

%build
if test "$USE_ICC"; then
%configure --enable-icc
else
%configure
fi
make %{?_smp_mflags}

%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%doc AUTHORS BUGS ChangeLog COPYRIGHT HISTORY INSTALL README THANKS doc/README.DOC doc/sextractor.pdf
%{_bindir}/sex
%{_bindir}/ldactoasc
%{_mandir}/man1/sex.1*
%{_mandir}/manx/sex.x*
%{_datadir}/sextractor

%changelog
* Thu Feb 25 2021 Emmanuel Bertin <bertin@iap.fr>
- Automatic RPM rebuild
* Tue May 13 2003 Emmanuel Bertin <bertin@iap.fr>
- RPM build for V2.3
* Fri Apr 04 2003 Emmanuel Bertin <bertin@iap.fr>
- RPM build for V2.3b4
* Wed Mar 05 2003 Emmanuel Bertin <bertin@iap.fr>
- RPM build for V2.3b3
* Fri Feb 07 2003 Emmanuel Bertin <bertin@iap.fr>
- Second RPM build
* Fri Jan 24 2003 Emmanuel Bertin <bertin@iap.fr>
- Second RPM build
* Sun Dec 15 2002 Emmanuel Bertin <bertin@iap.fr>
- First RPM build

# end of file
