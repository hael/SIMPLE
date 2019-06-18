!! Modified by Michael Eager, Feb 2018
program simple_test_img_export
include 'simple_lib.f08'
use simple_jpg, only : test_jpg_export
use simple_test_export_jpg, only : test_jpg_image
#ifdef USING_TIFF
use simple_tifflib
use simple_tifflib_test
#endif
use simple_cmdline, only: cmdline
implicit none
#include "simple_local_flags.inc"
type(cmdline)         :: cline
logical               :: be_verbose=.false.
character(len=8)      :: datestr
character(len=STDLEN) :: folder, command, olddir


if( command_argument_count() > 0 )then
   write(logfhandle,'(a)',advance='no') 'simple_test_img_export vol1=<volume.mrc> msk=<mask radius(in pixels)>'
   write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
   call cline%parse_oldschool
   call cline%check
   if( cline%defined('vol1') )then
      if( trim(cline%get_carg('vol1')) .eq. '' )then
         print *,' No input volume for testing'
      endif
   endif
   be_verbose = .false.
   if( cline%defined('verbose') )then
      if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
         be_verbose = .true.
      endif
   endif
endif

call seed_rnd
call date_and_time(date=datestr)
folder = trim('SIMPLE_TEST_IMGEXPORT_'//datestr)
call simple_mkdir(folder)
call simple_chdir(folder, olddir)

call create_raw_png_tmp
!call test_jpg_image(.true.)
!call test_jpg_export

#ifdef USING_TIFF
!call test_tiff_write1
!call test_tiff_write2
!call test_tiff_write3
!call test_tiff_write4

!call test_bigtiff_write
!call test_bigtiff_write1
!call test_bigtiff_write2

call test_bigtiff_write3

#endif

contains
        subroutine create_raw_png_tmp
        integer :: fid,ios
        integer(2) :: b
        integer(4) ::reclen

        fid=20
        inquire(iolength=reclen)b
        open(unit=20,file='bbb',form='unformatted',access='stream',status='replace',iostat=ios)
        !        open(NEWUNIT=fid,FILE="bbb.raw",ACCESS='STREAM',STATUS="REPLACE",IOSTAT=ios)!,FORM='UNFORMATTED')
        if (ios/=0) return
        write(fid) INT(Z'474e50890a1a0a0d0d00000052444849',16)
        write(fid) INT(Z'4b0000004b00000000000208ed2cb700',16)
        write(fid) INT(Z'000000bd4d416704b100004161fc0b8f',16)
        write(fid) INT(Z'00000005474b620600ff0044a0ff00ff',16)
        write(fid) INT(Z'0093a7bd700900000073594800120b00',16)
        write(fid) INT(Z'01120b00fc7eddd207000000454d4974',16)
        write(fid) INT(Z'040cd207952d38130075909049520f00',16)
        write(fid) INT(Z'785441445d5bad9c8db92393752904cc',16)
        write(fid) INT(Z'763d78cfe2f70fd8ecebffff67aedc70',16)
        write(fid) INT(Z'022aa5a7048240f7daef545942850ebe',16)
        write(fid) INT(Z'4aaaa51fa0489122e7f1f9aa1fd5e33f',16)
        write(fid) INT(Z'cf6b14a9a651a4a23414444f7f9abf59',16)
        write(fid) INT(Z'bda1ed35f587699597cd7f43a20d22c7',16)
        write(fid) INT(Z'f6e1dbb2c90622581c6f05bf181dbeb2',16)
        write(fid) INT(Z'a72ec036cd9660b5c45d3f0efd758e4f',16)
        write(fid) INT(Z'02c73715175c205f6445ec2977ec2af9',16)
        write(fid) INT(Z'5a72c08ee8e5d72230ec373ffa436c8b',16)
        write(fid) INT(Z'41db1162ee52d30dd5e15e10b721f99f',16)
        write(fid) INT(Z'2698deabb58e43bb0a346c6879f8c202',16)
        write(fid) INT(Z'ac788cc35bb238cf172389c4cb725784',16)
        write(fid) INT(Z'a5a6e576438cdbad08ba9f8a29dd8602',16)
        write(fid) INT(Z'fea21d86b486b684dc72eaff6dd6dc4d',16)
        write(fid) INT(Z'f7168426bf5a532288d92ecc86e48c47',16)
        write(fid) INT(Z'004b08e1b8b7dba9276df58563238983',16)
        write(fid) INT(Z'a8ea2287363269765de02fce5f55c581',16)
        write(fid) INT(Z'4f52a7e9156f123adc6df8cc7f514ce1',16)
        write(fid) INT(Z'b8ba99816cdb1bc03e4e3da3698d7877',16)
        write(fid) INT(Z'9ac428b88c9e016fb43ace434caed092',16)
        write(fid) INT(Z'783ca5932ea05faaeb145d21952b7525',16)
        write(fid) INT(Z'49d20cc06deebafca22a94b493e30432',16)
        write(fid) INT(Z'74523532bc2b88b3828d66456d046c0e',16)
        write(fid) INT(Z'351788bf75f241a200d75e2586e29931',16)
        write(fid) INT(Z'd52fe9a1516aec6845f1b8a3f550ba57',16)
        write(fid) INT(Z'858a1c5471eac79717c91c8b9e25db4d',16)
        write(fid) INT(Z'53c1a188d33d8617b1ccfde624e7eabd',16)
        write(fid) INT(Z'2aa851307a7e26ada890e9e2986d064b',16)
        write(fid) INT(Z'de02fcc83f3ad345edf53d656de05ca2',16)
        write(fid) INT(Z'4fd17238a2ea2a5984514a4650ef80bc',16)
        write(fid) INT(Z'3ae6eae437857c9c51511fd414c9e038',16)
        write(fid) INT(Z'b8ae11b7996c8bcd58e11ece13e4556c',16)
        write(fid) INT(Z'cd38e8f1db03387f532b2b1f521e22dc',16)
        write(fid) INT(Z'6475f45d4c06b6c5eb28cd1e4dbcfd0e',16)
        write(fid) INT(Z'1380a28275eb373334a4ff25de56e71d',16)
        write(fid) INT(Z'897cebfa481cd5dbc49719f3755c319b',16)
        write(fid) INT(Z'9933c6a5d1534c106e335883c548656e',16)
        write(fid) INT(Z'38cf8a1a9f57698222acc494fe9c188a',16)
        write(fid) INT(Z'a04894ddf19d5ca845450055a9ea1254',16)
        write(fid) INT(Z'9ea2cf48bb154d7b15ac02366add6acb',16)
        write(fid) INT(Z'170df9c1e3c51c6ec062455a8bd912f2',16)
        write(fid) INT(Z'e924f591bc2e766026493363d8be4902',16)
        write(fid) INT(Z'939f8d0920b576b10f06803c1a71d026',16)
        write(fid) INT(Z'e177b2698a348021b5cf076293a985f4',16)
        write(fid) INT(Z'de09cb410642a6cec99206bc8c51ee7a',16)
        write(fid) INT(Z'0ccc0902080357e704f1181235a7c0a8',16)
        write(fid) INT(Z'ee0b63b64c5cddc7de0afcb43f239ca4',16)
        write(fid) INT(Z'114bc15f85bf4b5e2392eea835e1cd3f',16)
        write(fid) INT(Z'40fe81a615e05b33aea3545e16e17dae',16)
        write(fid) INT(Z'ce115e3eaf3e4f754822627b64839e15',16)
        write(fid) INT(Z'4464c01a20d05fcfa0b660f0b47303b1',16)
        write(fid) INT(Z'ae4011e17c29a845f9341b1535960258',16)
        write(fid) INT(Z'b2a0caace17cfd04496294571c94b606',16)
        write(fid) INT(Z'fd78725b7354c105365ec30f0230b666',16)
        write(fid) INT(Z'd146812060048f065862a746f27ba9df',16)
        write(fid) INT(Z'df02627d1863a1ea5d37b1d6a35c3cb1',16)
        write(fid) INT(Z'0d99916f426a245e30477f6d997592db',16)
        write(fid) INT(Z'710ccce9002603831863009bd64d4ba0',16)
        write(fid) INT(Z'd24478b91d6e29e51d8c10862bd725c0',16)
        write(fid) INT(Z'332cfeaa395cfae0930a022c2a7a4440',16)
        write(fid) INT(Z'4ab12d564c011acb8b7c43d106d6f9fe',16)
        write(fid) INT(Z'53630868c1a55e9aa3f2f4bc77c64f76',16)
        write(fid) INT(Z'948fe781423360662de126119319ba59',16)
        write(fid) INT(Z'006d7a411044c2850035028619b0a291',16)
        write(fid) INT(Z'5e367b6ad148fa0b603ce56854ba7bb5',16)
        write(fid) INT(Z'a98b5dd524a9aca1d092a4a4d9823819',16)
        write(fid) INT(Z'a529066415b9367bcd18531389702aa8',16)
        write(fid) INT(Z'1331335764bb8c93d15a79ea84005bca',16)
        write(fid) INT(Z'0c5af04dcf7c160da20c1a55bf5b1832',16)
        write(fid) INT(Z'b5d34a1aba7d5c144b354f7721327283',16)
        write(fid) INT(Z'24d8adbe7ee12dc966aa4d578da895a2',16)
        write(fid) INT(Z'b8d6cdae99a99a8ad92ee321f47a71ca',16)
        write(fid) INT(Z'4ed21e0c6cd59851848e24bcac6978aa',16)
        write(fid) INT(Z'45e87f460c9c3662616dd57a1cc1e540',16)
        write(fid) INT(Z'6cb79514a18b4de9c8a9eb624424879a',16)
        write(fid) INT(Z'090b5d58154034925bd2de968d99a729',16)
        write(fid) INT(Z'd68a4588648b5d5c9c778598b4357e61',16)
        write(fid) INT(Z'5b4652f4e9f8b4f729c1660a53c13248',16)
        write(fid) INT(Z'199562432908824a46802326bb096b81',16)
        write(fid) INT(Z'6f1a9a72b0e56ee91b7c245aba1bb198',16)
        write(fid) INT(Z'417b7f8a10785dc10e1fda9b1a0b7696',16)
        write(fid) INT(Z'025a7152fab2a1e790c2f14fc96f124c',16)
        write(fid) INT(Z'b62ac1e111449024504e2f343e97aa08',16)
        write(fid) INT(Z'd6030aadc33353fc79cbdd2de1d0c6a2',16)
        write(fid) INT(Z'177af60d6e7571326c54a78bc0e3e34a',16)
        write(fid) INT(Z'59829a75563b2687aae577cbfb84a424',16)
        write(fid) INT(Z'd120169b29899092b08a2521dd0e1481',16)
        write(fid) INT(Z'4e0482a28513b98352dc84834c41672e',16)
        write(fid) INT(Z'65d5ab3fefac6a22a6d1872eb429eec6',16)
        write(fid) INT(Z'4e9ea3a12aca28f02788bd2d22092ef3',16)
        write(fid) INT(Z'a9aacf7296638ae99312449456e49250',16)
        write(fid) INT(Z'462c1515a06ac255337506603231430a',16)
        write(fid) INT(Z'0f599a67a03335a059f7b00da6c8de01',16)
        write(fid) INT(Z'eeb57cd06bbdd4698f3a160fb71ac310',16)
        write(fid) INT(Z'88a15a4e65169ba5f0e8f06a4b74a524',16)
        write(fid) INT(Z'fee92fb7adc922267e271f0ca613803c',16)
        write(fid) INT(Z'680d51d0c9649fa8c6c6050508d8eb4c',16)
        write(fid) INT(Z'08fad50520d552fbccc9a326a0a65743',16)
        write(fid) INT(Z'5ee7e9afeb6b7a5cc2809f8d8813525b',16)
        write(fid) INT(Z'ce22851ae52924cca7f5fcb713292924',16)
        write(fid) INT(Z'a18d48b3be437ee2a61f5fe33885140f',16)
        write(fid) INT(Z'50806117530ce40166656a33aad469c5',16)
        write(fid) INT(Z'9e5dad37118c88dd3978b74e45e19aa2',16)
        write(fid) INT(Z'ba9bd8138f6a27094df2c92bd436bd15',16)
        write(fid) INT(Z'd29493c429fd7d3de5b96531e52dfbf7',16)
        write(fid) INT(Z'ac2b251bf2cf3ea82fa7f43869e56138',16)
        write(fid) INT(Z'3751474504c119a3e8a0d9810dab5450',16)
        write(fid) INT(Z'7a866864a18bd84f52ccdea732b631b6',16)
        write(fid) INT(Z'e85e84f6e7eb4d2921888a182515d52b',16)
        write(fid) INT(Z'cc6aa45194924beddfdf96e45fef6fdf',16)
        write(fid) INT(Z'e617fbf3f3e7666cc987c7905f2cf1f3',16)
        write(fid) INT(Z'462cf1df66748cf40acc984202a3548c',16)
        write(fid) INT(Z'6c922f4a3b33769e4a1bd1efdb42b1e3',16)
        write(fid) INT(Z'cccf1a7a332e963f4b4907606d6ba735',16)
        write(fid) INT(Z'22915f6b49efdd39b7efde99b7d7fbdb',16)
        write(fid) INT(Z'7feedf6ff0fe96fdfcf9c7ae6ffc77e7',16)
        write(fid) INT(Z'1ca3862ad66e9a99008307f890819dd5',16)
        write(fid) INT(Z'97cd6a64c169a1bb6bb8d9ae6270860c',16)
        write(fid) INT(Z'a6677c54875b2c6c4b4c31e7a119d612',16)
        write(fid) INT(Z'71779b09142cc68ff94a5b91dcbfdf9e',16)
        write(fid) INT(Z'cbedfaff7f27fedf1e3ea3fb1cbdfdbc',16)
        write(fid) INT(Z'748e3cfa6522796a831061c3984614d1',16)
        write(fid) INT(Z'550840d23690bb4c439aa1c93ef0badd',16)
        write(fid) INT(Z'e96547f66026ea62170dafb3cfe9400f',16)
        write(fid) INT(Z'74a59194bebd62e339892053e7bc3fa7',16)
        write(fid) INT(Z'7fcc7f6ffafe93f99dfcff2396531bf7',16)
        write(fid) INT(Z'abbdea24462a19a9a9a4a66f8c09a2a8',16)
        write(fid) INT(Z'0bc70d4e05aa3e17a9bb9215b9eaa703',16)
        write(fid) INT(Z'77029b92a80ef070206ee9f66d402494',16)
        write(fid) INT(Z'f49e9af68f4c1b55fd0f1f5378fe87f5',16)
        write(fid) INT(Z'aa9a71e8c082dcc14cd52500008b8dd9',16)
        write(fid) INT(Z'f8979a5096ebeb5e2d13671bbcff3f5b',16)
        write(fid) INT(Z'63f627e830c41a0be661b866d2210029',16)
        write(fid) INT(Z'5ba481b48cfad5163cf9cb53cefcff9e',16)
        write(fid) INT(Z'cff3dfdb3f38f5dff1cbf1fe394b3cf1',16)
        write(fid) INT(Z'0a3c918b35c368148c494798c76bf156',16)
        write(fid) INT(Z'1db714d60c860d62734a1981b92b04f5',16)
        write(fid) INT(Z'c81a3a7ddb646f4942a4820040528248',16)
        write(fid) INT(Z'c3140580782ee0b30f1e9a9e1cbdfcfc',16)
        write(fid) INT(Z'd4cdfbcac7e78f8e2cf3c78f1c314147',16)
        write(fid) INT(Z'442c0505f4a8d40122da09168ba1b1d6',16)
        write(fid) INT(Z'983d7b2aa477785ccd4e1ebd2f353bbb',16)
        write(fid) INT(Z'68b181329b38ab11294042664a802490',16)
        write(fid) INT(Z'4b010195bf8cebed786f903e14790a1c',16)
        write(fid) INT(Z'3279f3c352ce0325e3f9e79e797e20a1',16)
        write(fid) INT(Z'409c54f0494c0a012210143543a358a5',16)
        write(fid) INT(Z'faa0de43ee1665ebd450d46dac9a7bb9',16)
        write(fid) INT(Z'8ff51f4cf4c20b1f5b0b7ad66533203b',16)
        write(fid) INT(Z'b6b6ac0638e21f6c7e07cbf053c782fa',16)
        write(fid) INT(Z'674cf28fc715086975fe20691389e3e0',16)
        write(fid) INT(Z'5027603c9083f4e03cda203857826d1a',16)
        write(fid) INT(Z'859d67b5a874ccbb2da70f2da5fb691a',16)
        write(fid) INT(Z'5b2d5eb9bbbbd20684025a4e66293782',16)
        write(fid) INT(Z'a4b50d04943500c1afb09d873b33705f',16)
        write(fid) INT(Z'8cc3b50f5fc05924c4f1f03e4e04f0d3',16)
        write(fid) INT(Z'402c22b098a44d018f06c0352e48a724',16)
        write(fid) INT(Z'9c1f414a7dd152b609c1b8a918b9bf63',16)
        write(fid) INT(Z'51c47ef6d7d5e88543344caa415f75aa',16)
        write(fid) INT(Z'81419a317b4f2862e59ba5aa70281464',16)
        write(fid) INT(Z'a7138e28e00f00e24b1a3c6950850059',16)
        write(fid) INT(Z'50d69c78c9d182f5c4664589483d3f22',16)
        write(fid) INT(Z'ddca6db5695ef4e0118b182e0679b027',16)
        write(fid) INT(Z'89b8628929259b6a62f4675006290628',16)
        write(fid) INT(Z'19b2d6ab869c42ec9ed13803901f580b',16)
        write(fid) INT(Z'1acf4c542b089f0463eb393a683efb52',16)
        write(fid) INT(Z'2f38ce540bc074be352a663e48e71523',16)
        write(fid) INT(Z'453265167c1817bae74cb13c8425b521',16)
        write(fid) INT(Z'04ecc592f0c6ba0b90d016340d9c7820',16)
        write(fid) INT(Z'b2d022adf3f59636c9eabc0b53e89d62',16)
        write(fid) INT(Z'b0abb58c97193c2b73cceb6c84e070c4',16)
        write(fid) INT(Z'86a4a416b5a29a9b1a50d503d5a766b4',16)
        write(fid) INT(Z'4e61568e3a275ac8f3c6411eb391b3b4',16)
        write(fid) INT(Z'55e1eaa14f3820a2af107b7ab307a1ae',16)
        write(fid) INT(Z'dd7a314db1cc3e93665d6bd65c8f0de1',16)
        write(fid) INT(Z'1840d1e056f4b4f553bf5a4130b5f6ae',16)
        write(fid) INT(Z'0f66d6837ab4152e197a3dea0f0ed3cf',16)
        write(fid) INT(Z'2e188d947cc60f0ff5a96d214b9e03a3',16)
        write(fid) INT(Z'c437a0a206893b9b628c3ce94ea18d7e',16)
        write(fid) INT(Z'e961ca8d62b7205e8db3466cadd43d46',16)
        write(fid) INT(Z'c79a8585c4b1b5105261c4f3fd97bebc',16)
        write(fid) INT(Z'978bc3a359aa03a5874b3a72dfabc128',16)
        write(fid) INT(Z'e553a8a2c82a2f58cf1bc9986f0b9203',16)
        write(fid) INT(Z'435243376f06f09810f447593544c493',16)
        write(fid) INT(Z'3d3c9c7b4b03354aa21524e375bc4280',16)
        write(fid) INT(Z'a6dce592c917ce307a9d68730de0d897',16)
        write(fid) INT(Z'26f5337f3879cc542f69485d950f254e',16)
        write(fid) INT(Z'240ce60c0ed68c830a02946d3cf66a05',16)
        write(fid) INT(Z'639ab4ed070b2d43f38791da58835b3c',16)
        write(fid) INT(Z'618e008486d981b1d879eba9c2ef71b3',16)
        write(fid) INT(Z'43004b9519974e2a22a8428842486281',16)
        write(fid) INT(Z'9cfcd6499eaa70290f192065a3e30c7a',16)
        write(fid) INT(Z'ea30cd5a4170ce5d71e167d88cd6bc22',16)
        write(fid) INT(Z'23ef0b8d86c04a74189fecf0d4b4688c',16)
        write(fid) INT(Z'b6534931b9c756e7422af1d6f8f64213',16)
        write(fid) INT(Z'866fc687964a2857602134ae2bc06bd7',16)
        write(fid) INT(Z'eb893ac28c6c7026cd816f76da290d5c',16)
        write(fid) INT(Z'c5ed4767b5fa7de8a89830873ca656f9',16)
        write(fid) INT(Z'86e88164e7d17d9fd094a40c40041b10',16)
        write(fid) INT(Z'51b092be32fe37b4b7efa311aa40b36d',16)
        write(fid) INT(Z'6056ab776983d9ed8c2588f39ce7a3b1',16)
        write(fid) INT(Z'bfa2f610bf56d787e119facc0febc801',16)
        write(fid) INT(Z'4b96f47af9b2f34e60b4c9dbf67d0dc4',16)
        write(fid) INT(Z'fa1923196c086a74eaf3145716da6b5f',16)
        write(fid) INT(Z'7d09ef5eafd2d6d69769ed44758466e5',16)
        write(fid) INT(Z'aa0c2c3d21671f16aa51619231842382',16)
        write(fid) INT(Z'14b073f0f536563544368eb231d273af',16)
        write(fid) INT(Z'e039b24f78365a45f008a28db6834b2a',16)
        write(fid) INT(Z'50e4a25eb87c2e375307b7394e8f4de2',16)
        write(fid) INT(Z'a756c74832b31b9dd395da289b3acb7a',16)
        write(fid) INT(Z'12cb4e9003c0ba28e6d3dde0612e04a6',16)
        write(fid) INT(Z'7cf85d9c908306cbbf011e6cc1b60977',16)
        write(fid) INT(Z'844a89369f27b9775e1c54434a1bc1d9',16)
        write(fid) INT(Z'a17babb3fba844782988037f12b3343e',16)
        write(fid) INT(Z'd330dc150db1e2c026d1f4c1296284c6',16)
        write(fid) INT(Z'8617360c237a3af29378170f0eab854d',16)
        write(fid) INT(Z'dc5c025a298a45dced13f13ade198a1e',16)
        write(fid) INT(Z'7ee764682e4e19afc74ccd37fd55e19a',16)
        write(fid) INT(Z'28bbc3df49dd38cdea821457e426fd71',16)
        write(fid) INT(Z'5db1e0cf51a8badc4c972844d5c60d5b',16)
        write(fid) INT(Z'0f75ad5125ef0dbf575366c2e4eb1737',16)
        write(fid) INT(Z'fad5cdd6d5c53d792f67ae404d77966d',16)
        write(fid) INT(Z'fc21dd132d66f0fd776e858db6859437',16)
        write(fid) INT(Z'c49ce3a545bc2b003c15ed7e7db6e60a',16)
        write(fid) INT(Z'e0253b9b9fdeae2ce879ebc3dc5390c5',16)
        write(fid) INT(Z'aabccc963b78086bc3cc5c1b18e91c67',16)
        write(fid) INT(Z'6fb87630af0bfa88dfb106ec8610d704',16)
        write(fid) INT(Z'175784eb2fa43d7ef84852f688b1bc39',16)
        write(fid) INT(Z'14fbd73e31cca8de3baa0368b085dfb2',16)
        write(fid) INT(Z'ea902dc5f71f5cd35ae30e73b6299afb',16)
        write(fid) INT(Z'1b13877d9df4e097c01eac2ce886a644',16)
        write(fid) INT(Z'b98787df8f14bc5a70ddbdeb6c611d8d',16)
        write(fid) INT(Z'bc49790af546665de3146a8a996f0b56',16)
        write(fid) INT(Z'd3f0f52cd5c216ff006aecca3414f763',16)
        write(fid) INT(Z'06761bc6bea80eea1e7b8b0ad6f9ac9f',16)
        write(fid) INT(Z'fe96a785d0eff315adc20d9578f1372e',16)
        write(fid) INT(Z'97961ff93f1319a60f8ea2a3edd9659b',16)
        write(fid) INT(Z'40a64247fe7a5a56f0f007ff5c5a7baa',16)
        write(fid) INT(Z'a6f0e7f804bbf2a43dd7b63c82d6022c',16)
        write(fid) INT(Z'f9d1e0b33c31ffe8e97d2cacd96ca7c2',16)
        write(fid) INT(Z'3b74eba49a44684ecbaf597116dd65c6',16)
        write(fid) INT(Z'ff0ab886d9003c0d7c8635247b76affa',16)
        write(fid) INT(Z'b2378a3562e3965ccede5b56e4b1fdff',16)
        write(fid) INT(Z'ff88bae1e18f3c255f3a9ff6d35ed4d8',16)
        write(fid) INT(Z'17e2c3ae7ac9fae9fff6655dff801e19',16)
        write(fid) INT(Z'8a5bea039e9c3b830000009e4e454900',16)
        write(fid) INT(Z'6042ae44',4), INT(Z'0082',2)
        close(fid)
    end subroutine create_raw_png_tmp

end program simple_test_img_export
